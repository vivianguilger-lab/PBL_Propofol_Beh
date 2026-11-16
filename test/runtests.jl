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
