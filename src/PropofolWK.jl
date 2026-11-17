# Módulo principal: junta os arquivos do modelo e diz quais funções ficam
# disponíveis para quem fizer `using PropofolWK`.
module PropofolWK

# LinearAlgebra: contas com matrizes (usada na farmacocinética)
using LinearAlgebra
# Random: sorteios reprodutíveis (usada nos pacientes virtuais)
using Random

include("farmacocinetica.jl")
include("alvo.jl")
include("hemodinamica.jl")
include("risco.jl")
include("populacao.jl")
include("grafico.jl")

export parametros_eleveld, simular_pk, taxa_bolus, taxa_alvo
export ParametrosSu, simular_hemodinamica
export LIMIARES, nivel, percentil, faixa_pam, periodos_de_risco, descrever
export ω²_ELEVELD, ω²_SU_C50_RPT, sortear_paciente, simular_populacao
export grafico_faixa

end
