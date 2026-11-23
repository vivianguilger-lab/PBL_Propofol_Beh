# =============================================================================
# servidor.jl — a plataforma: um pequeno servidor web local
#
# O médico abre a página no navegador (app/index.html), preenche os dados e
# clica em "Prever". A página pede a previsão a este servidor, que roda no
# próprio computador (endereço 127.0.0.1: nada sai da máquina), calcula com
# `prever` e devolve o gráfico e os períodos de risco.
#
# Feito só com o pacote Sockets, que já vem com o Julia: não precisa instalar
# nada da internet. Atende uma pessoa por vez, o que basta para uso local.
#
# Endereços que o servidor entende:
#   /                  → a página (app/index.html)
#   /api/prever?...    → a previsão, em JSON (formato de dados que a página lê)
#   /api/comparar?...  → pressão × dose de propofol (remifentanil fixo): paciente típico, 61 doses ou alvos
#   /api/remi?...      → pressão × dose de remifentanil (propofol fixo): paciente típico, 26 taxas
# =============================================================================

# ---------------------------------------------------------------- JSON
# Transforma valores do Julia em texto JSON. Só os tipos que a resposta usa.
json(x::AbstractString) = "\"" * escape_json(x) * "\""
json(x::Bool) = x ? "true" : "false"
json(::Nothing) = "null"
json(x::Symbol) = json(string(x))
json(x::Integer) = string(x)
# Números com 4 algarismos significativos: 98,49 mmHg nas curvas, mas 0,0327 nos parâmetros (com 2 casas virava 0,03).
json(x::Real) = isfinite(x) ? string(round(x; sigdigits = 4)) : "null"
json(x::AbstractVector) = "[" * join(json.(x), ",") * "]"
json(x::NamedTuple) = "{" * join(["\"$k\":" * json(v) for (k, v) in pairs(x)], ",") * "}"

# No JSON, aspas, barras invertidas e quebras de linha precisam de "\" antes.
escape_json(s) = replace(s, "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n", "\r" => "\\r", "\t" => "\\t")

# ---------------------------------------------------------------- endereço
# Desfaz a codificação de endereço da web: "+" vira espaço e "%2C" vira ",".
function decodificar_url(s)
    c = codeunits(replace(s, "+" => " "))     # os bytes do texto
    bytes = UInt8[]
    i = 1
    while i <= length(c)
        if c[i] == UInt8('%') && i + 2 <= length(c)
            push!(bytes, parse(UInt8, String(c[i+1:i+2]); base = 16))   # "%2C" → byte 0x2C
            i += 3
        else
            push!(bytes, c[i])
            i += 1
        end
    end
    return String(bytes)
end

# "idade=70&peso=80" → Dict("idade" => "70", "peso" => "80")
function ler_parametros(consulta)
    p = Dict{String, String}()
    for par in split(consulta, "&"; keepempty = false)
        chave, valor = occursin("=", par) ? split(par, "="; limit = 2) : (par, "")
        p[decodificar_url(chave)] = decodificar_url(valor)
    end
    return p
end

# ---------------------------------------------------------------- previsão
# Lê um número do formulário; aceita vírgula decimal ("2,5"). Erro claro se faltar.
function numero(p, chave, nome)
    haskey(p, chave) && !isempty(strip(p[chave])) || error("Preencha o campo: $nome.")
    v = tryparse(Float64, replace(strip(p[chave]), "," => "."))
    v === nothing && error("Valor inválido em $nome: \"$(p[chave])\".")
    return v
end

# Lê os dados do paciente, a entrada do propofol (dose ou alvo) e o remifentanil dos campos do formulário.
# A manutenção chega em µg/kg/min (como na bomba) e vai para o modelo em mg/kg/h (× 60 ÷ 1000).
# O remifentanil chega em µg/kg/min e o início dele em minutos antes do propofol (vai em segundos).
function ler_formulario(p)
    dc = haskey(p, "dc") && !isempty(strip(p["dc"])) ? numero(p, "dc", "débito cardíaco") : nothing
    paciente = (idade = numero(p, "idade", "idade"), peso = numero(p, "peso", "peso"),
                altura = numero(p, "altura", "altura"), sexo = get(p, "sexo", ""),
                pas = numero(p, "pas", "PAS"), pad = numero(p, "pad", "PAD"),
                pam = numero(p, "pam", "PAM"), fc = numero(p, "fc", "FC"), dc = dc)
    modo = get(p, "modo", "dose")
    manut = haskey(p, "manut") && !isempty(strip(p["manut"])) ? numero(p, "manut", "manutenção") : 0.0
    entrada = modo == "alvo" ? (alvo_ugml = numero(p, "alvo", "concentração-alvo"),) :
        (dose_mgkg = numero(p, "dose", "dose"),
         duracao_bolus_s = round(Int, numero(p, "bolus", "tempo de injeção")),
         manutencao_mgkgh = manut * 60 / 1000)
    remi = haskey(p, "remi") && !isempty(strip(p["remi"])) ? numero(p, "remi", "remifentanil") : 0.0
    antes = haskey(p, "remi_antes") && !isempty(strip(p["remi_antes"])) ? numero(p, "remi_antes", "início do remifentanil") : 9.0
    return paciente, modo, entrada, (remi_ugkgmin = remi, remi_antes_s = round(Int, antes * 60))
end

# Resposta JSON com tratamento de erro: se `f()` der erro, devolve {"ok": false, "erro": ...}.
function responder(f)
    try
        return json(f())
    catch e
        msg = e isa ErrorException ? e.msg : "Não foi possível calcular: " * sprint(showerror, e)
        return json((ok = false, erro = msg))
    end
end

# Tira uma amostra a cada `passo` segundos (a página não precisa dos 1201 pontos).
amostrar(v; passo = 2) = v[1:passo:end]

"""
    responder_previsao(p)

Recebe os campos do formulário (Dict) e devolve o texto JSON da previsão: as
curvas (tempo, faixas da PAM, da FC e da Cp, RPT e VS), o resumo (`metricas`),
os períodos de risco, a dose total e a farmacocinética e os parâmetros usados;
ou `{"ok": false, "erro": "mensagem"}`.
"""
responder_previsao(p) = responder() do
    paciente, modo, entrada, remi = ler_formulario(p)
    r = prever(; paciente..., entrada..., remi...)
    a = amostrar
    curvas = (t = a(r.t), pam_inf = a(r.faixa.inf), pam_med = a(r.faixa.mediana), pam_sup = a(r.faixa.sup),
              fc_inf = a(r.fc.inf), fc_med = a(r.fc.mediana), fc_sup = a(r.fc.sup),
              cp_inf = a(r.cp.inf), cp_med = a(r.cp.mediana), cp_sup = a(r.cp.sup),
              rpt = a(100 .* (r.rpt .- 1)), vs = a(100 .* (r.vs .- 1)),
              fcp = a(100 .* (r.fc.mediana ./ paciente.fc .- 1)),
              pam_q1 = a(r.metade.inf), pam_q3 = a(r.metade.sup),
              bis_inf = a(r.bis.inf), bis_q1 = a(r.metade_bis.inf), bis_med = a(r.bis.mediana),
              bis_q3 = a(r.metade_bis.sup), bis_sup = a(r.bis.sup), dormindo = a(r.dormindo),
              cr_inf = a(r.cr.inf), cr_med = a(r.cr.mediana), cr_sup = a(r.cr.sup),
              pas_inf = a(r.pas.inf), pas_med = a(r.pas.mediana), pas_sup = a(r.pas.sup),
              pad_inf = a(r.pad.inf), pad_med = a(r.pad.mediana), pad_sup = a(r.pad.sup),
              dc_inf = a(r.dc.inf), dc_med = a(r.dc.mediana), dc_sup = a(r.dc.sup))
    # Windkessel do paciente e onda de pressão do paciente típico (a cada 5 ms)
    wk = (Rp = r.wk.Rp, C = r.wk.C, Rc = r.wk.Rc, L = r.wk.L, dc = r.wk.dc, vs = r.wk.vs, tau = round(r.wk.Rp * r.wk.C; digits = 2),
          pas = r.wk.pas, pad = r.wk.pad, pam = r.wk.pam, pas_braco = paciente.pas)
    ondas = [(t = o.t, pressao = o.pressao, fluxo = o.fluxo, pas = o.pas, pad = o.pad, pam = o.pam, fc = o.fc) for o in r.ondas]
    periodos = [(inicio = x.inicio, fim = x.fim, pior = x.pior, melhor = x.melhor,
                 pam_min = x.pam_min, pam_max = x.pam_max) for x in r.periodos]
    su = ParametrosSu()
    parametros = (c50_rpt = su.c50_rpt, emax_rpt = su.emax_rpt, γ_rpt = su.γ_rpt, c50_vs = su.c50_vs,
                  emax_vs = PropofolWK.emax_vs(su, paciente.idade), k_out = su.k_out * 60, fb = su.fb,
                  ω²_c50 = ω²_SU_C50_RPT, n = 500, c50_remi_rpt = su.c50_remi_rpt, emax_remi_rpt = su.emax_remi_rpt,
                  sl_remi_fc = su.sl_remi_fc, sl_remi_vs = su.sl_remi_vs, int_rpt = su.int_rpt, int_vs = su.int_vs,
                  int_fc = su.int_fc, ω²_emax_remi = ω²_SU_EMAX_REMI_RPT, α = su.α, c50_int_fc = su.c50_int_fc,
                  ω²_eleveld = ω²_ELEVELD, ω²_minto = ω²_MINTO, opioide = remi.remi_ugkgmin > 0)
    bis = ParametrosBIS()
    sono = (r.sono..., ce50 = PropofolWK.ce50_bis(bis, paciente.idade), ke0 = PropofolWK.ke0_bis(bis, paciente.peso) * 60,
            atraso = PropofolWK.atraso_bis(bis, paciente.idade), ω²_ce50 = ω²_ELEVELD_BIS.ce50, ω²_ke0 = ω²_ELEVELD_BIS.ke0)
    (ok = true, modo = modo, dose_mg = sum(r.taxa), remi_ug = r.taxa_remi === nothing ? 0.0 : sum(r.taxa_remi),
     remi = remi, curvas = curvas, metricas = r.metricas, sono = sono, wk = wk, ondas = ondas,
     periodos = periodos, texto = r.texto, pk = r.pk, pk_remi = r.pk_remi, parametros = parametros)
end

"""
    responder_comparar(p)

Gráfico pressão × dose de propofol, com o remifentanil fixo no valor do formulário:
o mesmo paciente (paciente típico) com 61 doses, de 0 a 3,0 mg/kg de 0,05 em 0,05,
ou 61 concentrações-alvo, de 0 a 6 µg/ml de 0,1 em 0,1, conforme o modo. Para cada uma: PAM mínima, queda, segundos abaixo de 65 e
de 55, o nível de risco e o BIS mínimo (ver `risco_por_dose`).
"""
responder_comparar(p) = responder() do
    paciente, modo, entrada, remi = ler_formulario(p)
    valores = modo == "alvo" ? collect(0.0:0.1:6.0) : collect(0.0:0.05:3.0)
    extra = modo == "alvo" ? (;) : (duracao_bolus_s = entrada.duracao_bolus_s, manutencao_mgkgh = entrada.manutencao_mgkgh)
    linhas = risco_por_dose(valores; paciente..., modo, extra..., remi...)
    (ok = true, modo = modo, linhas = linhas)
end

"""
    responder_remi(p)

Gráfico pressão × dose de remifentanil, com o propofol fixo no valor do formulário
(dose ou alvo): o paciente típico com 26 taxas de remifentanil, de 0 a 0,5
µg/kg/min de 0,02 em 0,02 (a maior taxa média nas 21 induções do VitalDB foi
0,46). Mesmas colunas de `responder_comparar` (ver `risco_por_remi`).
"""
responder_remi(p) = responder() do
    paciente, modo, entrada, remi = ler_formulario(p)
    linhas = risco_por_remi(collect(0.0:0.02:0.5); paciente..., entrada..., remi_antes_s = remi.remi_antes_s)
    (ok = true, modo = modo, linhas = linhas)
end

# ---------------------------------------------------------------- servidor
# Página da plataforma (fica em app/index.html, na raiz do projeto).
const PAGINA = joinpath(@__DIR__, "..", "app", "index.html")

# Escreve uma resposta HTTP completa (cabeçalho + conteúdo) na conexão.
function enviar(sock, status, tipo, corpo)
    write(sock, "HTTP/1.1 $status\r\nContent-Type: $tipo; charset=utf-8\r\n" *
                "Content-Length: $(sizeof(corpo))\r\nConnection: close\r\n\r\n", corpo)
end

# Atende um pedido do navegador: lê a primeira linha ("GET /caminho HTTP/1.1"),
# ignora o resto do cabeçalho e responde conforme o caminho.
function atender(sock)
    linha = readline(sock)
    while !isempty(readline(sock)) end          # pula o resto do cabeçalho
    partes = split(linha)
    alvo = length(partes) >= 2 ? partes[2] : "/"
    caminho, consulta = occursin("?", alvo) ? split(alvo, "?"; limit = 2) : (alvo, "")
    if caminho == "/"
        enviar(sock, "200 OK", "text/html", read(PAGINA, String))
    elseif caminho == "/api/prever"
        enviar(sock, "200 OK", "application/json", responder_previsao(ler_parametros(consulta)))
    elseif caminho == "/api/comparar"
        enviar(sock, "200 OK", "application/json", responder_comparar(ler_parametros(consulta)))
    elseif caminho == "/api/remi"
        enviar(sock, "200 OK", "application/json", responder_remi(ler_parametros(consulta)))
    else
        enviar(sock, "404 Not Found", "text/plain", "Página não encontrada.")
    end
    close(sock)
end

"""
    iniciar_servidor(porta)

Começa a atender em http://127.0.0.1:porta, em segundo plano. Devolve o
servidor; `close(servidor)` para.
"""
function iniciar_servidor(porta)
    servidor = listen(IPv4(127, 0, 0, 1), porta)
    @async while isopen(servidor)
        sock = try accept(servidor) catch; break end
        try
            atender(sock)
        catch e
            @warn "Erro ao atender um pedido" exception = e
            close(sock)
        end
    end
    return servidor
end

# Abre o endereço no navegador padrão (Windows, macOS ou Linux).
function abrir_navegador(url)
    try
        Sys.iswindows() ? run(`cmd /c start $url`) :
        Sys.isapple() ? run(`open $url`) : run(`xdg-open $url`)
    catch
        println("Abra no navegador: $url")
    end
end

"""
    plataforma(; porta = 8080, abrir = true)

Abre a plataforma: inicia o servidor e o navegador. Fica rodando até a janela
do Julia ser fechada (ou Ctrl+C).
"""
function plataforma(; porta = 8080, abrir = true)
    servidor = iniciar_servidor(porta)
    url = "http://127.0.0.1:$porta"
    println("Plataforma aberta em $url  (para fechar: feche esta janela ou aperte Ctrl+C)")
    abrir && abrir_navegador(url)
    try
        while isopen(servidor)
            sleep(1)
        end
    catch e
        e isa InterruptException || rethrow()
    finally
        close(servidor)
    end
end
