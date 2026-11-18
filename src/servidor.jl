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
# =============================================================================

# ---------------------------------------------------------------- JSON
# Transforma valores do Julia em texto JSON. Só os tipos que a resposta usa.
json(x::AbstractString) = "\"" * escape_json(x) * "\""
json(x::Bool) = x ? "true" : "false"
json(x::Symbol) = json(string(x))
json(x::Integer) = string(x)
json(x::Real) = isfinite(x) ? string(round(x; digits = 2)) : "null"
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

"""
    responder_previsao(p)

Recebe os campos do formulário (Dict) e devolve o texto JSON da resposta:
`{"ok": true, "svg": ..., "texto": [...], "periodos": [...]}` ou
`{"ok": false, "erro": "mensagem"}`.
"""
function responder_previsao(p)
    try
        paciente = (idade = numero(p, "idade", "idade"), peso = numero(p, "peso", "peso"),
                    altura = numero(p, "altura", "altura"), sexo = get(p, "sexo", ""),
                    pam = numero(p, "pam", "PAM"), fc = numero(p, "fc", "FC"))
        modo = get(p, "modo", "dose")
        entrada = modo == "alvo" ? (alvo_ugml = numero(p, "alvo", "concentração-alvo"),) :
            (dose_mgkg = numero(p, "dose", "dose"),
             duracao_bolus_s = round(Int, numero(p, "bolus", "tempo de injeção")),
             manutencao_mgkgh = haskey(p, "manut") && !isempty(strip(p["manut"])) ? numero(p, "manut", "manutenção") : 0.0)
        r = prever(; paciente..., entrada...)
        titulo = modo == "alvo" ? "Alvo plasmático de $(p["alvo"]) µg/ml" :
                 "$(p["dose"]) mg/kg em $(p["bolus"]) s"
        periodos = [(inicio = x.inicio, fim = x.fim, pior = x.pior, melhor = x.melhor,
                     pam_min = x.pam_min, pam_max = x.pam_max) for x in r.periodos]
        return json((ok = true, svg = svg_faixa(r.t, r.faixa; titulo), texto = r.texto, periodos = periodos))
    catch e
        msg = e isa ErrorException ? e.msg : "Não foi possível calcular: " * sprint(showerror, e)
        return json((ok = false, erro = msg))
    end
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
