# =============================================================================
# exemplos/prever_paciente.jl — exemplo de uso: previsão para um paciente
#
# Rodar no Julia, com o projeto ativado (] activate .):
#     include("exemplos/prever_paciente.jl")
# Os gráficos (.svg) são gravados na pasta resultados/; abra no navegador.
# =============================================================================
using PropofolWK

# Dados do paciente: troque aqui para simular outro paciente.
paciente = (idade = 70, peso = 80, altura = 175, sexo = "M", pas = 100.0, pad = 62.0, pam = 75.0, fc = 70.0)

# As duas entradas que o médico pode escolher: nome do arquivo, entrada e título.
entradas = [
    ("dose", (dose_mgkg = 2.0, duracao_bolus_s = 30), "2 mg/kg em 30 s"),
    ("alvo", (alvo_ugml = 3.0,), "alvo plasmático de 3 µg/ml"),
]

# Pasta resultados/ na raiz do projeto (criada se não existir).
pasta = joinpath(@__DIR__, "..", "resultados")
mkpath(pasta)

# Descrição do paciente para o título do gráfico.
quem = "$(paciente.sexo == "M" ? "Homem" : "Mulher"), $(paciente.idade) anos, PA basal $(round(Int, paciente.pas))/$(round(Int, paciente.pad)) mmHg (PAM $(round(Int, paciente.pam)))"

for (nome, entrada, titulo) in entradas
    r = prever(; paciente..., entrada...)
    arquivo = grafico_faixa(joinpath(pasta, "faixa_$nome.svg"), r.t, r.faixa; titulo = "$quem — $titulo")
    println("\n== $titulo")
    foreach(println, r.texto)
    println("gráfico: ", normpath(arquivo))
end
