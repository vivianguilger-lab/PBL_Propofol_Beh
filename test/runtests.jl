# Testes do PropofolWK. Rodar no Julia, na pasta do projeto:  ] test
using PropofolWK
using Test, Random, Sockets

# ler_csv(arquivo): lê uma tabela gerada pelo PAS. A 1ª linha tem os nomes das
# colunas; as outras, números separados por vírgula. Devolve uma função col:
# col("cp") é a coluna "cp" inteira.
function ler_csv(arquivo)
    linhas = readlines(joinpath(@__DIR__, arquivo))
    cab = split(linhas[1], ",")
    valores = [parse.(Float64, split(l, ",")) for l in linhas[2:end]]   # uma lista de números por linha
    return nome -> [v[findfirst(==(nome), cab)] for v in valores]
end

# Tabela gerada com test/gerar_referencia_pk.py: Cp do modelo de Eleveld no PAS,
# para dois pacientes (A e B), com e sem opioide.
@testset "Eleveld 2018: farmacocinética" begin
    col = ler_csv("referencia_pk.csv")
    pacientes = (A = (62, 72, 160, "F"), B = (30, 80, 180, "M"))
    for (nome, (idade, peso, altura, sexo)) in pairs(pacientes)
        for opioide in (false, true)
            pk = parametros_eleveld(idade, peso, altura, sexo; opioide)
            cp = simular_pk(pk, col("taxa_$nome"))
            ref = col("cp_$(nome)_$(opioide ? "com" : "sem")_opioide")
            # diferença menor que 0,0001 µg/ml em todos os 1200 segundos
            @test maximum(abs.(cp .- ref)) < 1e-4
        end
    end
end

# Tabela gerada com test/gerar_referencia_pas.py: a Cp de entrada (bolus de
# 1,5 mg/kg em 30 s, PK de Eleveld do PAS, paciente de 62 anos) e a PAM e a FC
# que o PAS calcula.
@testset "Su 2023, só propofol" begin
    col = ler_csv("referencia_pas.csv")
    r = simular_hemodinamica(col("cp"), 104.0, 78.0, 62)

    @testset "antes do propofol, tudo no basal" begin
        @test r.pam[1] ≈ 104.0
        @test r.fc[1] ≈ 78.0
    end
    @testset "igual ao PAS (diferença < 0,001 mmHg e 0,001 bpm)" begin
        @test maximum(abs.(r.pam .- col("pam_pas"))) < 1e-3
        @test maximum(abs.(r.fc .- col("fc_pas"))) < 1e-3
        @test maximum(abs.(r.vs_rel .- col("vs_rel_exato"))) < 1e-5
        @test maximum(abs.(r.rpt_rel .- col("rpt_rel_exato"))) < 1e-5
    end
    @testset "PAM não depende do VS basal" begin
        # no PAS, VS basal de 70 e de 90 ml deram a mesma PAM; aqui o VS basal nem entra
        @test minimum(r.pam) ≈ minimum(col("pam_pas")) atol = 1e-3
    end
end

@testset "Faixa e períodos de risco" begin
    @testset "níveis de risco" begin
        @test nivel(70.0) == :baixo
        @test nivel(65.0) == :baixo       # 65 já é baixo (PAM ≥ 65)
        @test nivel(60.0) == :moderado
        @test nivel(55.0) == :moderado
        @test nivel(50.0) == :alto
    end
    @testset "percentil (igual ao PERCENTIL do Excel)" begin
        # no Excel: =PERCENTIL({10;20;30;40};0,25) dá 17,5
        @test percentil([40.0, 10.0, 30.0, 20.0], 0.25) ≈ 17.5
        @test percentil([40.0, 10.0, 30.0, 20.0], 0.5) ≈ 25.0
        @test percentil([40.0, 10.0, 30.0, 20.0], 1.0) ≈ 40.0
    end
    @testset "faixa: percentis entre curvas" begin
        # 101 curvas constantes, de 0 a 100 mmHg: percentil 5 = 5, mediana = 50, percentil 95 = 95
        curvas = [fill(Float64(v), 4) for v in 0:100]
        f = faixa_pam(curvas)
        @test f.inf ≈ fill(5.0, 4)
        @test f.mediana ≈ fill(50.0, 4)
        @test f.sup ≈ fill(95.0, 4)
    end
    @testset "períodos" begin
        # faixa inventada só para testar a lógica: 6 instantes, de 0 a 5 s
        t = 0:5
        inf = [80.0, 70.0, 60.0, 50.0, 60.0, 70.0]
        sup = [90.0, 85.0, 75.0, 70.0, 75.0, 80.0]
        p = periodos_de_risco(t, inf, sup)
        @test [(x.inicio, x.fim, x.pior, x.melhor) for x in p] ==
              [(0, 1, :baixo, :baixo), (2, 2, :moderado, :baixo), (3, 3, :alto, :baixo),
               (4, 4, :moderado, :baixo), (5, 5, :baixo, :baixo)]
        @test p[1].pam_min == 70.0 && p[1].pam_max == 90.0
        @test descrever(p)[3] == "de 3 a 3 s: PAM entre 50 e 70 mmHg — risco entre alto e baixo"
    end
end

@testset "Pacientes virtuais (variabilidade publicada)" begin
    taxa = taxa_bolus(1.5 * 72, 30, 301)       # 300 s, para o teste ser rápido
    paciente = (62, 72, 160, "F", 104.0, 78.0)
    wk = parametrizar_windkessel(pam = 104.0, pad = 82.0, fc = 78.0, peso = 72, altura = 160)
    @testset "mesma semente, mesmas curvas" begin
        a = simular_populacao(paciente..., taxa; wk, n = 5, semente = 7)
        b = simular_populacao(paciente..., taxa; wk, n = 5, semente = 7)
        c = simular_populacao(paciente..., taxa; wk, n = 5, semente = 8)
        @test a.curvas == b.curvas
        @test a.curvas != c.curvas
    end
    @testset "sem variabilidade, todas as versões são o paciente típico" begin
        zero_pk = map(_ -> 0.0, ω²_ELEVELD)    # todos os ω² = 0
        pop = simular_populacao(paciente..., taxa; wk, n = 3, ω²_pk = zero_pk, ω²_c50 = 0.0, ω²_bis = (ce50 = 0.0, ke0 = 0.0))
        tipico = simular_circulacao(simular_pk(parametros_eleveld(62, 72, 160, "F"), taxa), nothing, 62, wk;
                                    pam0 = 104.0, fc0 = 78.0).pam
        @test all(c -> c ≈ tipico, pop.curvas)
        bis_tipico = simular_bis(parametros_eleveld(62, 72, 160, "F"), taxa, 62, 72).bis
        @test all(b -> b ≈ bis_tipico, pop.bis)
    end
    @testset "o sorteio tem média 0 e variância ω² (em escala log)" begin
        rng = Xoshiro(1)
        pk = parametros_eleveld(62, 72, 160, "F")
        psu = ParametrosSu()
        versoes = [sortear_paciente(pk, psu, rng) for _ in 1:20000]
        η_v1 = [log(v[1].V1 / pk.V1) for v in versoes]
        η_c50 = [log(v[2].c50_rpt / psu.c50_rpt) for v in versoes]
        media(x) = sum(x) / length(x)
        variancia(x) = sum((x .- media(x)) .^ 2) / (length(x) - 1)
        @test abs(media(η_v1)) < 0.02
        @test variancia(η_v1) ≈ ω²_ELEVELD.V1 rtol = 0.05
        @test abs(media(η_c50)) < 0.02
        @test variancia(η_c50) ≈ ω²_SU_C50_RPT rtol = 0.05
    end
end

@testset "Arredondamento evita períodos picotados" begin
    # borda de cima "raspando" os 65 mmHg: 65,2 / 64,8 / 65,3 / 64,9
    t = 0:3
    inf = [60.0, 60.0, 60.0, 60.0]
    sup = [65.2, 64.8, 65.3, 64.9]
    @test length(periodos_de_risco(t, inf, sup; arredondar = false)) == 4   # picotado
    @test length(periodos_de_risco(t, inf, sup)) == 1                       # tudo vira 65: um período só
end

@testset "Gráfico da faixa (SVG)" begin
    @test PropofolWK.escala(0, 0, 10, 100, 200) == 100      # início do intervalo
    @test PropofolWK.escala(10, 0, 10, 100, 200) == 200     # fim do intervalo
    @test PropofolWK.escala(5, 0, 10, 100, 200) == 150      # meio
    t = collect(0.0:600.0)
    faixa = (inf = fill(60.0, 601), mediana = fill(70.0, 601), sup = fill(80.0, 601))
    arquivo = grafico_faixa(joinpath(mktempdir(), "teste.svg"), t, faixa; titulo = "Teste <&>")
    svg = read(arquivo, String)
    @test startswith(svg, "<svg")
    @test occursin("Tempo após o início do bolus (s)", svg)
    @test occursin("PAM (mmHg)", svg)
    @test occursin("<polygon", svg) && occursin("<polyline", svg)
    @test occursin("Teste &lt;&amp;&gt;", svg)              # título com caracteres especiais escapados
end

@testset "Concentração-alvo (TCI no plasma)" begin
    pk = parametros_eleveld(62, 72, 160, "F")
    @testset "sem limite de bomba: Cp igual ao alvo a partir do 2º segundo" begin
        taxa = taxa_alvo(pk, 3.0, 600)
        cp = simular_pk(pk, taxa)
        @test cp[1] == 0.0                      # antes da infusão
        @test all(c -> c ≈ 3.0, cp[2:end])      # depois, sempre no alvo
        @test all(>=(0), taxa)                  # nunca negativa
    end
    @testset "com limite de bomba: sobe até o alvo e fica nele" begin
        taxa = taxa_alvo(pk, 3.0, 600; taxa_max = 2.0)
        cp = simular_pk(pk, taxa)
        @test maximum(taxa) <= 2.0
        @test cp[5] < 3.0                       # nos primeiros segundos ainda não chegou
        @test cp[end] ≈ 3.0
    end
    @testset "alvo que cai: a bomba para até a Cp baixar" begin
        alvos = [fill(4.0, 300); fill(2.0, 300)]
        taxa = taxa_alvo(pk, alvos, 600)
        cp = simular_pk(pk, taxa)
        @test taxa[301] == 0.0                  # logo após a troca, bomba parada
        @test cp[end] ≈ 2.0                     # e depois mantém o novo alvo
    end
end

@testset "Previsão completa" begin
    paciente = (idade = 70, peso = 80, altura = 175, sexo = "M", pas = 105.0, pad = 60.0, pam = 75.0, fc = 70.0)
    @testset "modo dose: igual a montar as partes à mão" begin
        r = prever(; paciente..., dose_mgkg = 2.0, duracao_s = 300, n = 20)
        wk = parametrizar_windkessel(pam = 75.0, pad = 60.0, fc = 70.0, peso = 80, altura = 175)
        pop = simular_populacao(70, 80, 175, "M", 75.0, 70.0, taxa_bolus(160.0, 30, 301); wk, n = 20)
        @test r.faixa.mediana ≈ faixa_pam(pop.curvas).mediana
        @test sum(r.taxa) ≈ 160.0               # 2 mg/kg × 80 kg
    end
    @testset "modo dose com manutenção" begin
        r = prever(; paciente..., dose_mgkg = 2.0, manutencao_mgkgh = 6.0, duracao_s = 300, n = 5)
        @test r.taxa[1] ≈ 160.0 / 30            # durante o bolus
        @test r.taxa[end] ≈ 6.0 * 80 / 3600     # depois, só a manutenção (mg/s)
    end
    @testset "modo concentração-alvo" begin
        r = prever(; paciente..., alvo_ugml = 3.0, duracao_s = 300, n = 5)
        @test r.taxa ≈ taxa_alvo(parametros_eleveld(70, 80, 175, "M"), 3.0, 301)
    end
    @testset "entradas inválidas dão erro com explicação" begin
        @test_throws ErrorException prever(; paciente...)                                   # nenhuma entrada
        @test_throws ErrorException prever(; paciente..., dose_mgkg = 2.0, alvo_ugml = 3.0) # as duas
        @test_throws ErrorException prever(; paciente..., sexo = "X", dose_mgkg = 2.0)
        @test_throws ErrorException prever(; paciente..., peso = -1, dose_mgkg = 2.0)
        @test_throws ErrorException prever(; paciente..., pad = 80.0, dose_mgkg = 2.0)       # PAD acima da PAM
    end
end

@testset "Plataforma (servidor local)" begin
    @testset "leitura do endereço" begin
        @test PropofolWK.decodificar_url("2%2C5+mg") == "2,5 mg"
        @test PropofolWK.ler_parametros("idade=70&sexo=M&manut=") == Dict("idade" => "70", "sexo" => "M", "manut" => "")
    end
    @testset "JSON" begin
        @test PropofolWK.json((ok = true, x = 1.23456, s = "a\"b", v = [1, 2])) == "{\"ok\":true,\"x\":1.235,\"s\":\"a\\\"b\",\"v\":[1,2]}"
        @test PropofolWK.json(0.0327) == "0.0327"           # parâmetros pequenos não perdem algarismos
    end
    # pede um endereço ao servidor, como o navegador faz, e devolve o conteúdo da resposta
    function pedir(porta, caminho)
        sock = connect(IPv4(127, 0, 0, 1), porta)
        write(sock, "GET $caminho HTTP/1.1\r\nHost: teste\r\n\r\n")
        resposta = read(sock, String)
        close(sock)
        return split(resposta, "\r\n\r\n"; limit = 2)
    end
    servidor = PropofolWK.iniciar_servidor(8765)
    try
        cab, corpo = pedir(8765, "/")
        @test startswith(cab, "HTTP/1.1 200")
        @test occursin("<form", corpo)
        cab, corpo = pedir(8765, "/api/prever?idade=70&sexo=M&peso=80&altura=175&pas=105&pad=60&pam=75&fc=70&modo=dose&dose=2&bolus=30&manut=")
        @test occursin("\"ok\":true", corpo)
        @test occursin("\"pam_med\":[", corpo)                # curvas para os gráficos da página
        @test occursin("\"risco\":", corpo)                   # resumo do cartão de risco
        @test occursin("\"dormindo\":[", corpo)              # % dormindo, para a faixa do sono
        @test occursin("\"sono\":{", corpo)                  # quando, em média, dorme e acorda
        cab, corpo = pedir(8765, "/api/prever?idade=70&sexo=M&peso=80&altura=175&pas=105&pad=60&pam=75&fc=70&modo=alvo&alvo=3%2C0")
        @test occursin("\"ok\":true", corpo)                  # vírgula decimal aceita
        cab, corpo = pedir(8765, "/api/prever?idade=70&sexo=M&peso=&altura=175&pas=105&pad=60&pam=75&fc=70&modo=dose&dose=2&bolus=30")
        @test occursin("Preencha o campo: peso.", corpo)     # erro explicado para o médico
        cab, corpo = pedir(8765, "/api/comparar?idade=70&sexo=M&peso=80&altura=175&pas=105&pad=60&pam=75&fc=70&modo=dose&dose=2&bolus=30&manut=0")
        @test occursin("\"linhas\":[", corpo)                 # gráfico pressão × dose (cartão 04)
        @test count("\"valor\":", corpo) == 61                # 0 a 3,0 mg/kg, de 0,05 em 0,05
        cab, corpo = pedir(8765, "/api/remi?idade=70&sexo=M&peso=80&altura=175&pas=105&pad=60&pam=75&fc=70&modo=dose&dose=2&bolus=30&manut=0&remi=0.1&remi_antes=9")
        @test count("\"valor\":", corpo) == 26                # remifentanil de 0 a 0,5 µg/kg/min, de 0,02 em 0,02
        cab, corpo = pedir(8765, "/api/prever?idade=70&sexo=M&peso=80&altura=175&pas=105&pad=60&pam=75&fc=70&modo=dose&dose=2&bolus=30&manut=0&remi=0.1&remi_antes=9")
        @test occursin("\"cr_med\":[", corpo)                # Cp do remifentanil, para o gráfico
        cab, corpo = pedir(8765, "/outra")
        @test startswith(cab, "HTTP/1.1 404")
    finally
        close(servidor)
    end
end

@testset "Risco por dose (paciente típico)" begin
    paciente = (idade = 81, peso = 64, altura = 168, sexo = "M", pas = 168.0, pad = 76.0, pam = 107.0, fc = 64.0)
    r = risco_por_dose([1.0, 2.0, 3.0]; paciente...)
    @test risco_por_dose([0.0]; paciente...)[1].pam_min ≈ 107.0 atol = 1e-6   # sem propofol, a PAM fica na basal
    # igual a simular o paciente típico à mão, dose por dose
    pk = parametros_eleveld(81, 64, 168, "M")
    wk = parametrizar_windkessel(pam = 107.0, pad = 76.0, fc = 64.0, peso = 64, altura = 168)
    manual = minimum(simular_circulacao(simular_pk(pk, taxa_bolus(2.0 * 64, 30, 1201)), nothing, 81, wk;
                                        pam0 = 107.0, fc0 = 64.0).pam)
    @test r[2].pam_min ≈ manual
    @test r[1].pam_min > r[2].pam_min > r[3].pam_min    # mais dose, PAM mínima mais baixa
    @test r[2].risco == nivel(round(manual))
    @test r[1].bis_min > r[2].bis_min > r[3].bis_min    # mais dose, BIS mais baixo
end

# Tabela gerada com test/gerar_referencia_bis.py: Ce e BIS de Eleveld no PAS, sem opioide,
# para os mesmos dois pacientes da farmacocinética.
@testset "Eleveld 2018: BIS" begin
    col = ler_csv("referencia_bis.csv")
    pacientes = (A = (62, 72, 160, "F"), B = (30, 80, 180, "M"))
    for (nome, (idade, peso, altura, sexo)) in pairs(pacientes)
        r = simular_bis(parametros_eleveld(idade, peso, altura, sexo), col("taxa_$nome"), idade, peso)
        @test maximum(abs.(r.ce .- col("ce_$nome"))) < 1e-6     # concentração no cérebro igual ao PAS
        # a equação do BIS é igual à do PAS: com o mesmo atraso que o PAS usa (inteiro, uma amostra
        # a menos), a diferença é zero
        p = ParametrosBIS()
        n = ceil(Int, PropofolWK.atraso_bis(p, idade)) - 1
        calc = [PropofolWK.bis_de_ce(c, PropofolWK.ce50_bis(p, idade), p) for c in r.ce]
        @test maximum(abs.([k > n ? calc[k - n] : 93.0 for k in eachindex(calc)] .- col("bis_$nome"))) < 1e-6
        # com o atraso exato do artigo, a diferença para o PAS fica abaixo de 1 ponto de BIS
        @test maximum(abs.(r.bis .- col("bis_$nome"))) < 1.0
    end
    @testset "acordado e limiar" begin
        p = ParametrosBIS()
        @test PropofolWK.bis_de_ce(0.0, 3.0, p) == 93.0          # sem propofol, BIS de acordado
        @test PropofolWK.bis_de_ce(3.0, 3.0, p) ≈ 46.5           # no Ce50, metade do efeito máximo
        m = metricas_bis(collect(0.0:4.0), [93.0, 70.0, 55.0, 50.0, 65.0])
        @test m.bis_min == 50.0 && m.t_hipnose == 2.0 && m.s_hipnose == 2.0 && m.acordou
    end
end

# Tabela gerada com test/gerar_referencia_remi.py: Cp do remifentanil (Minto), Cp do propofol
# (Eleveld com opioide) e PAM e FC do Su com os dois fármacos, no PAS. Remifentanil de 0,12
# µg/kg/min desde 540 s antes do propofol; propofol 1,5 mg/kg em 30 s; mulher de 62 anos.
@testset "Remifentanil (Minto 1997 e Su 2023)" begin
    col = ler_csv("referencia_remi.csv")
    pk = parametros_eleveld(62, 72, 160, "F"; opioide = true)
    pr = parametros_minto(62, 72, 160, "F")
    tp = vcat(zeros(540), taxa_bolus(1.5 * 72, 30, 1201))
    tr = taxa_remifentanil(72; ugkgmin = 0.12, antes_s = 540, duracao_s = 1200)
    cp, cr = simular_pk(pk, tp), simular_pk(pr, tr)
    @test maximum(abs.(cr .- col("cr"))) < 1e-6              # Cp do remifentanil igual à do PAS (ng/ml)
    @test maximum(abs.(cp .- col("cp"))) < 1e-6              # propofol com opioide igual ao PAS
    h = simular_hemodinamica(cp, 104.0, 78.0, 62; cr)
    @test maximum(abs.(h.pam .- col("pam"))) < 1e-3          # PAM igual à do PAS
    @test maximum(abs.(h.fc .- col("fc"))) < 1e-3            # FC igual à do PAS
    # sem remifentanil, igual à versão só com propofol
    cp0 = simular_pk(parametros_eleveld(62, 72, 160, "F"), taxa_bolus(1.5 * 72, 30, 1201))
    @test simular_hemodinamica(cp0, 104.0, 78.0, 62; cr = zeros(1201)).pam == simular_hemodinamica(cp0, 104.0, 78.0, 62).pam
    @testset "pressão por dose de remifentanil" begin
        paciente = (idade = 62, peso = 72, altura = 160, sexo = "F", pas = 148.0, pad = 82.0, pam = 104.0, fc = 78.0)
        r = risco_por_remi([0.0, 0.1, 0.2]; paciente..., dose_mgkg = 1.5)
        @test r[1].pam_min > r[2].pam_min > r[3].pam_min    # mais remifentanil, PAM mínima mais baixa
        @test r[1].pam_min ≈ risco_por_dose([1.5]; paciente...)[1].pam_min   # remi 0 = só propofol
        @test risco_por_dose([1.5]; paciente..., remi_ugkgmin = 0.1)[1].pam_min ≈ r[2].pam_min
    end
    @testset "previsão com remifentanil" begin
        r = prever(idade = 62, peso = 72, altura = 160, sexo = "F", pas = 148.0, pad = 82.0, pam = 104.0, fc = 78.0, dose_mgkg = 1.5,
                   remi_ugkgmin = 0.12, n = 50)
        @test r.faixa.mediana[1] < 104.0                     # no início do propofol, o remi já baixou a PAM
        @test r.cr.mediana[1] > 0                            # e já há remifentanil no sangue
        @test_throws ErrorException prever(idade = 62, peso = 72, altura = 160, sexo = "F", pas = 148.0, pad = 82.0, pam = 104.0, fc = 78.0,
                                           dose_mgkg = 1.5, remi_ugkgmin = -0.1)
    end
    @testset "sorteio do remifentanil: média 0 e covariância publicadas" begin
        rng = Xoshiro(2)
        psu, prm = ParametrosSu(), parametros_minto(62, 72, 160, "F")
        v = [sortear_remifentanil(prm, psu, rng) for _ in 1:20000]
        e = [x[2].emax_remi_rpt - psu.emax_remi_rpt for x in v]
        a = [x[2].sl_remi_fc - psu.sl_remi_fc for x in v]
        b = [x[2].sl_remi_vs - psu.sl_remi_vs for x in v]
        m(x) = sum(x) / length(x)
        cov(x, y) = sum((x .- m(x)) .* (y .- m(y))) / (length(x) - 1)
        @test abs(m(e)) < 0.02
        @test cov(e, e) ≈ ω²_SU_EMAX_REMI_RPT rtol = 0.05
        @test cov(a, a) ≈ 0.00382 rtol = 0.05
        @test cov(b, b) ≈ 0.00868 rtol = 0.05
        @test cov(a, b) ≈ 0.00329 rtol = 0.05
        @test cov(log.([x[1].CL for x in v]), log.([x[1].CL for x in v])) ≈ ω²_MINTO.CL rtol = 0.05
    end
end

# O centro do modelo: Windkessel de 4 elementos (Stergiopulos 1999), parametrizado com o
# monitor e acoplado ao modelo de Su em malha fechada.
@testset "Windkessel de 4 elementos" begin
    wk = parametrizar_windkessel(pam = 104.0, pad = 82.0, fc = 78.0, peso = 72, altura = 160)
    @testset "sem fluxo, a pressão cai exponencialmente: P(t) = P0·e^(−t/RpC)" begin
        Rp, C = 1.0, 1.2
        b = batimento(Rp, C, 0.056, 0.0051, 60.0, 0.0, (C * 100.0, 0.0); guardar = true)   # 1 s sem ejeção
        t = 0.005 .* (0:length(b.pressao) - 1)
        @test maximum(abs.(b.pressao .- 100.0 .* exp.(-t ./ (Rp * C)))) < 1e-9
    end
    @testset "PAM = Rp × fluxo médio (lei de Ohm da circulação)" begin
        @test wk.pam ≈ wk.Rp * wk.dc * 1000 / 60 rtol = 1e-9
        @test wk.pam ≈ 104.0 rtol = 1e-9                     # Rp = PAM/DC reproduz a PAM medida
        @test wk.pad ≈ 82.0 atol = 1e-6                       # C identificada pela PAD medida
        @test wk.pas > wk.pam > wk.pad
        @test 0.5 < wk.C < 2.0                                # perto da de Stergiopulos (1,22 mL/mmHg)
    end
    @testset "batimento periódico: o estado se repete" begin
        b = batimento(wk.Rp, wk.C, wk.Rc, wk.L, 78.0, wk.vs, wk.x0)
        @test all(abs.(b.x .- wk.x0) .< 1e-9)
    end
    @testset "malha fechada sem fármaco: tudo fica no basal" begin
        s = simular_circulacao(zeros(300), nothing, 62, wk; pam0 = 104.0, fc0 = 78.0)
        @test maximum(abs.(s.pam .- 104.0)) < 1e-9
        @test maximum(abs.(s.pad .- 82.0)) < 1e-6
    end
    @testset "congruência: a PAM do Windkessel acompanha a do modelo de Su" begin
        pk = parametros_eleveld(62, 72, 160, "F"; opioide = true)
        tp = vcat(zeros(540), taxa_bolus(1.5 * 72, 30, 1201))
        cr = simular_pk(parametros_minto(62, 72, 160, "F"), taxa_remifentanil(72; ugkgmin = 0.12, antes_s = 540, duracao_s = 1200))
        cp = simular_pk(pk, tp)
        su = simular_hemodinamica(cp, 104.0, 78.0, 62; cr).pam
        s = simular_circulacao(cp, cr, 62, wk; pam0 = 104.0, fc0 = 78.0)
        @test maximum(abs.(s.pam .- su)) < 1.0                # diferença só nos transitórios (< 1 mmHg)
        @test abs(minimum(s.pam) - minimum(su)) < 0.2
        @test all(s.pas .> s.pam .> s.pad)
        @test s.dc[end] < s.dc[1]                             # o débito cai com os fármacos
    end
end
