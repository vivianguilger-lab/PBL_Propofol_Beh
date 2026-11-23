# =============================================================================
# bis.jl — profundidade anestésica (BIS) durante a indução com propofol
#
# Modelo: Eleveld DJ, Colin P, Absalom AR, Struys MMRF. Pharmacokinetic–
# pharmacodynamic model for propofol for broad application in anaesthesia and
# sedation. Br J Anaesth 2018;120(5):942–59. doi:10.1016/j.bja.2018.01.018
# (o mesmo artigo da farmacocinética de farmacocinetica.jl).
#
# Implementação de referência: Python Anesthesia Simulator (PAS), commit 4afc03a:
#   src/python_anesthesia_simulator/pd_models.py, classe BIS_model, L205–L236
#   (parâmetros), L367–L392 (cálculo e atraso do monitor);
#   src/python_anesthesia_simulator/pk_models.py, L335 (ke0, sangue arterial).
#
# Só propofol: com remifentanil = 0, o modelo do PAS se reduz ao de Eleveld.
# Unidades: tempo em segundos; concentrações em µg/ml; BIS de 0 a 100.
# =============================================================================

# Parâmetros publicados por Eleveld 2018, como estão no PAS.
Base.@kwdef struct ParametrosBIS
    ce50_35::Float64     = 3.08      # Ce com metade do efeito aos 35 anos (µg/ml)       pd_models L222
    k_idade::Float64     = -0.00635  # mudança do Ce50 por ano de idade (exponencial)    pd_models L222
    γ_abaixo::Float64    = 1.89      # inclinação da curva com Ce ≤ Ce50                 pd_models L226
    γ_acima::Float64     = 1.47      # inclinação da curva com Ce > Ce50                 pd_models L228
    e0::Float64          = 93.0      # BIS acordado (sem propofol)                       pd_models L232
    ke0_70::Float64      = 0.146     # passagem sangue → cérebro, 70 kg (1/min)          pk_models L335
    atraso_base::Float64 = 15.0      # atraso do monitor: 15 + exp(k_atraso × idade) s   pd_models L219
    k_atraso::Float64    = 0.0517    #                                                   pd_models L236
end

# Limiar de hipnose adequada: BIS abaixo de 60. Definição clínica: a faixa-alvo do
# BIS na anestesia geral é de 40 a 60 (NICE, Depth of anaesthesia monitors, DG6, seção 4).
const LIMIAR_BIS = 60.0

# Ce50 do paciente: cai com a idade (PAS L218 e L222).
ce50_bis(p::ParametrosBIS, idade) = p.ce50_35 * exp(p.k_idade * (idade - 35))
# ke0 do paciente (1/s): menor em quem pesa mais (PAS pk_models L335; ÷ 60 para 1/s).
ke0_bis(p::ParametrosBIS, peso) = p.ke0_70 * (peso / 70)^(-0.25) / 60
# Atraso do monitor de BIS (s): maior em quem é mais velho (PAS L219 e L236).
atraso_bis(p::ParametrosBIS, idade) = p.atraso_base + exp(p.k_atraso * idade)

# BIS para uma concentração no cérebro: E0 − Emax × U^γ / (1 + U^γ), com U = Ce/Ce50
# e Emax = E0 (PAS L234 e L363). A inclinação γ muda acima do Ce50 (PAS L341).
function bis_de_ce(ce, ce50, p::ParametrosBIS)
    γ = ce <= ce50 ? p.γ_abaixo : p.γ_acima
    u = (ce / ce50)^γ
    return p.e0 - p.e0 * u / (1 + u)
end

"""
    simular_ce(pk, ke0, taxa; dt = 1.0)

Concentração no sítio efetor (cérebro), em µg/ml, a cada `dt` segundos. É um 4º
compartimento que recebe propofol do sangue com velocidade `ke0` (1/s) e não muda
a concentração do sangue, como no PAS (pk_models.py L533–L537).
"""
function simular_ce(pk, ke0, taxa; dt = 1.0)
    # taxas de transferência (1/min, divididas por 60 abaixo); a 4ª linha é o sítio efetor
    k10, k12, k13 = pk.CL / pk.V1, pk.Q2 / pk.V1, pk.Q3 / pk.V1
    k21, k31 = pk.Q2 / pk.V2, pk.Q3 / pk.V3
    A = [-(k10 + k12 + k13) / 60 k12 / 60 k13 / 60 0.0;
         k21 / 60 -k21 / 60 0.0 0.0;
         k31 / 60 0.0 -k31 / 60 0.0;
         ke0 0.0 0.0 -ke0]
    B = [1 / pk.V1, 0.0, 0.0, 0.0]
    Φ = exp(A .* dt)
    Γ = A \ ((Φ - I) * B)
    x = zeros(4)
    ce = zeros(length(taxa))
    for k in eachindex(taxa)
        ce[k] = x[4]                               # Ce no início do segundo k
        x = Φ * x .+ Γ .* taxa[k]
    end
    return ce
end

"""
    simular_bis(pk, taxa, idade, peso; p = ParametrosBIS(), dt = 1.0)

BIS como o monitor mostra, a cada `dt` segundos, para a infusão `taxa` (mg/s).
Devolve `(ce, bis)`. O monitor mostra com atraso (Eleveld 2018): o BIS do instante t
é o calculado no instante t − atraso, com o atraso exato do artigo (interpolando
entre duas amostras); antes do primeiro valor atrasado, o BIS de acordado.
Obs.: o PAS arredonda esse atraso para um número inteiro de amostras e usa uma a
menos (pd_models.py L309 e L385–L388); a diferença é menor que 1 s.
"""
function simular_bis(pk, taxa, idade, peso; p = ParametrosBIS(), dt = 1.0)
    ce = simular_ce(pk, ke0_bis(p, peso), taxa; dt)
    ce50 = ce50_bis(p, idade)
    calc = [bis_de_ce(c, ce50, p) for c in ce]        # BIS sem o atraso
    τ = atraso_bis(p, idade) / dt                      # atraso, em amostras (pode ser fracionário)
    bis = map(eachindex(calc)) do k
        s = (k - 1) - τ                                # posição atrasada (0 = primeira amostra)
        s < 0 && return p.e0
        i = floor(Int, s)
        f = s - i                                      # fração entre a amostra i e a seguinte
        i + 2 > length(calc) ? calc[i + 1] : (1 - f) * calc[i + 1] + f * calc[i + 2]
    end
    return (ce = ce, bis = bis)
end

"""
    metricas_bis(t, bis; limiar = LIMIAR_BIS)

Resumo da hipnose: BIS mínimo, instante em que o BIS fica abaixo do limiar pela
primeira vez (`t_hipnose`, s; `nothing` se nunca fica), segundos abaixo do limiar
(`s_hipnose`) e se no fim da simulação ele já voltou para cima do limiar (`acordou`).
"""
function metricas_bis(t, bis; limiar = LIMIAR_BIS)
    i = findfirst(<(limiar), bis)
    dt = t[2] - t[1]
    return (bis_min = minimum(bis), t_hipnose = i === nothing ? nothing : t[i],
            s_hipnose = dt * count(<(limiar), bis), acordou = i !== nothing && bis[end] >= limiar)
end
