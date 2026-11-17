# Testes do PropofolWK. Rodar no Julia, na pasta do projeto:  ] test
using PropofolWK
using Test

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
