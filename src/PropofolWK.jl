# Módulo principal: junta os arquivos do modelo e diz quais funções ficam
# disponíveis para quem fizer `using PropofolWK`.
module PropofolWK

# LinearAlgebra: contas com matrizes (usada na farmacocinética)
using LinearAlgebra

include("farmacocinetica.jl")
include("hemodinamica.jl")
include("risco.jl")

export parametros_eleveld, simular_pk, taxa_bolus
export ParametrosSu, simular_hemodinamica
export LIMIARES, nivel, percentil, faixa_pam, periodos_de_risco, descrever

end
