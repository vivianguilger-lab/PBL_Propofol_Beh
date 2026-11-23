# Módulo principal: junta os arquivos do modelo e diz quais funções ficam
# disponíveis para quem fizer `using PropofolWK`.
module PropofolWK

# LinearAlgebra: contas com matrizes (usada na farmacocinética)
using LinearAlgebra
# Random: sorteios reprodutíveis (usada nos pacientes virtuais)
using Random
# Sockets: conexões de rede (usada no servidor da plataforma)
using Sockets

include("farmacocinetica.jl")
include("alvo.jl")
include("remifentanil.jl")
include("bis.jl")
include("hemodinamica.jl")
include("risco.jl")
include("populacao.jl")
include("grafico.jl")
include("previsao.jl")
include("servidor.jl")

export parametros_eleveld, simular_pk, taxa_bolus, taxa_alvo
export parametros_minto, ω²_MINTO, taxa_remifentanil
export ParametrosBIS, LIMIAR_BIS, simular_bis, metricas_bis
export ParametrosSu, simular_hemodinamica
export LIMIARES, nivel, percentil, faixa_pam, periodos_de_risco, descrever
export ω²_ELEVELD, ω²_SU_C50_RPT, ω²_ELEVELD_BIS, ω²_SU_EMAX_REMI_RPT, sortear_paciente, sortear_remifentanil, simular_populacao
export svg_faixa, grafico_faixa
export prever, risco_por_dose, risco_por_remi, plataforma

end
