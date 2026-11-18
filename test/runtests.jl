# Testes do PropofolWK. Rodar no Julia, na pasta do projeto:  ] test
using PropofolWK
using Test, Random

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
    @testset "mesma semente, mesmas curvas" begin
        a = simular_populacao(paciente..., taxa; n = 5, semente = 7)
        b = simular_populacao(paciente..., taxa; n = 5, semente = 7)
        c = simular_populacao(paciente..., taxa; n = 5, semente = 8)
        @test a.curvas == b.curvas
        @test a.curvas != c.curvas
    end
    @testset "sem variabilidade, todas as versões são o paciente típico" begin
        zero_pk = map(_ -> 0.0, ω²_ELEVELD)    # todos os ω² = 0
        pop = simular_populacao(paciente..., taxa; n = 3, ω²_pk = zero_pk, ω²_c50 = 0.0)
        tipico = simular_hemodinamica(simular_pk(parametros_eleveld(62, 72, 160, "F"), taxa), 104.0, 78.0, 62).pam
        @test all(c -> c ≈ tipico, pop.curvas)
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
    paciente = (idade = 70, peso = 80, altura = 175, sexo = "M", pam = 75.0, fc = 70.0)
    @testset "modo dose: igual a montar as partes à mão" begin
        r = prever(; paciente..., dose_mgkg = 2.0, duracao_s = 300, n = 20)
        pop = simular_populacao(70, 80, 175, "M", 75.0, 70.0, taxa_bolus(160.0, 30, 301); n = 20)
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
    end
end
