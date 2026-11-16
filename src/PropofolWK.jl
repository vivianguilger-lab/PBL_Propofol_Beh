# Módulo principal: junta os arquivos do modelo e diz quais funções ficam
# disponíveis para quem fizer `using PropofolWK`.
module PropofolWK

include("hemodinamica.jl")

export ParametrosSu, simular_hemodinamica

end
