# Testes do PropofolWK. Rodar no Julia, na pasta do projeto:  ] test
using PropofolWK
using Test

# Tabela gerada rodando o próprio PAS (commit 4afc03a) com
# test/gerar_referencia_pas.py: a Cp de entrada (bolus de 1,5 mg/kg em
# 30 s, PK de Eleveld do PAS, paciente de 62 anos) e a PAM e a FC que o PAS calcula.
# Lê o CSV: a 1ª linha tem os nomes das colunas; as outras, números separados por vírgula.
linhas = readlines(joinpath(@__DIR__, "referencia_pas.csv"))
cab = split(linhas[1], ",")
valores = [parse.(Float64, split(l, ",")) for l in linhas[2:end]]   # uma lista de números por linha
# col("cp") devolve a coluna "cp" inteira
col(nome) = [v[findfirst(==(nome), cab)] for v in valores]
cp = col("cp")

@testset "Su 2023, só propofol" begin
    r = simular_hemodinamica(cp, 104.0, 78.0, 62)

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
