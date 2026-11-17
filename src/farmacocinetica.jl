# =============================================================================
# farmacocinetica.jl — da dose de propofol à concentração no sangue (Cp)
#
# Modelo: Eleveld DJ, Colin P, Absalom AR, Struys MMRF. Pharmacokinetic–
# pharmacodynamic model for propofol for broad application in anaesthesia and
# sedation. Br J Anaesth 2018;120(5):942–59. doi:10.1016/j.bja.2018.01.018
#
# Implementação de referência: Python Anesthesia Simulator (PAS), commit 4afc03a,
# src/python_anesthesia_simulator/pk_models.py, L256–L335 (modelo 'Eleveld').
#
# Escolhas: concentração arterial (o padrão do PAS; o artigo também tem uma
# variante para sangue venoso, não usada aqui); sem opioide por padrão, porque
# o projeto é só propofol.
#
# Unidades: tempo em segundos; dose em mg; volumes em L; Cp em mg/L = µg/ml.
# =============================================================================

# Valores θ do modelo de Eleveld, com a numeração do artigo (PAS L269–L287).
# Ficam de fora o θ7 (erro residual do ajuste) e o θ17 e θ18 (sangue venoso).
#                        dose (mg/s)
#                            ↓
#               músculo ⇄ SANGUE ⇄ gordura
#                 (V2)      (V1)      (V3)
#                            ↓
#                    eliminação (CL, fígado)
const θ_ELEVELD = (
    V1 = 6.2830780766822,     # θ1  V1 de referência (L)                       L270            Aqui permanece a mesma dúvida de se esses valores
    V2 = 25.5013145036879,    # θ2  V2 de referência (L)                       L271            deveriam ficar ou ser trocados 
    V3 = 272.8166615043603,   # θ3  V3 de referência (L)                       L272
    CL_m = 1.7895836588902,   # θ4  depuração de referência, homem (L/min)     L273
    Q2 = 1.7500983738779,     # θ5  Q2 de referência (L/min)                   L274
    Q3 = 1.1085424008536,     # θ6  Q3 de referência (L/min)                   L275
    CL_mat50 = 42.2760190602615,  # θ8  maturação da depuração: idade pós-menstrual com 50% (semanas)  L277
    CL_matγ = 9.0548452392807,    # θ9  inclinação dessa maturação                                     L278
    V2_idade = -0.015633,     # θ10 V2 menor com a idade (por ano)             L279
    CL_opioide = -0.00285709, # θ11 depuração menor com a idade, com opioide   L280
    V1_peso50 = 33.5531248778544, # θ12 peso com 50% do V1 máximo (kg)         L281
    V3_opioide = -0.0138166,  # θ13 V3 menor com a idade, com opioide          L282
    Q3_mat50 = 68.2767978846832,  # θ14 maturação de Q3 (semanas)              L283
    CL_f = 2.1002218877899,   # θ15 depuração de referência, mulher (L/min)    L284
    Q2_mat = 1.3042680471360, # θ16 Q2 maior enquanto Q3 amadurece             L285
)

# Paciente de referência do modelo: homem, 35 anos, 70 kg, 170 cm (PAS L262–L267).
const REF_ELEVELD = (idade = 35.0, peso = 70.0, altura = 170.0)

# Curva sigmoide usada pelo modelo (PAS L291).
sig(x, c50, γ) = x^γ / (c50^γ + x^γ)

# Massa livre de gordura de Al-Sallami, que o modelo usa para o V3 (PAS L294–L298).
function massa_magra_alsallami(masculino, peso, idade, imc)
    masculino ?
        (0.88 + (1 - 0.88) / (1 + (idade / 13.4)^(-12.7))) * (9270 * peso) / (6680 + 216 * imc) :
        (1.11 + (1 - 1.11) / (1 + (idade / 7.1)^(-1.1))) * (9270 * peso) / (8780 + 244 * imc)
end

"""
    parametros_eleveld(idade, peso, altura, sexo; opioide = false)

Volumes (L) e depurações (L/min) do modelo de Eleveld para um paciente.
`sexo` é "M" ou "F"; `altura` em cm. Com `opioide = true`, aplica os termos que
o artigo usa quando há opioide junto (depuração e V3 menores com a idade).
"""
function parametros_eleveld(idade, peso, altura, sexo; opioide = false)
    θ, ref = θ_ELEVELD, REF_ELEVELD
    masculino = sexo == "M"
    # idade pós-menstrual em anos (nascido a termo: + 40 semanas) e IMC (PAS L300–L301)
    ipm, ipm_ref = idade + 40 / 52, ref.idade + 40 / 52
    imc, imc_ref = peso / (altura / 100)^2, ref.peso / (ref.altura / 100)^2
    # maturação (em adultos ≈ 1, mas fica para seguir o artigo) (PAS L303–L306)
    mat_cl = sig(ipm * 52, θ.CL_mat50, θ.CL_matγ) / sig(ipm_ref * 52, θ.CL_mat50, θ.CL_matγ)
    mat_q3 = sig(ipm * 52, θ.Q3_mat50, 1) / sig(ipm_ref * 52, θ.Q3_mat50, 1)
    q3_imaturo = sig(ipm * 52, θ.Q3_mat50, 1)
    # massa magra do paciente / da referência (PAS L307–L308)
    mm = massa_magra_alsallami(masculino, peso, idade, imc) / massa_magra_alsallami(true, ref.peso, ref.idade, imc_ref)
    # efeito do opioide: exp(θ × idade) se houver, 1 se não houver (PAS L310–L313)
    opi(x) = opioide ? exp(x * idade) : 1.0
    # volumes (PAS L317, L320, L322)
    V1 = θ.V1 * sig(peso, θ.V1_peso50, 1) / sig(ref.peso, θ.V1_peso50, 1)
    V2 = θ.V2 * peso / ref.peso * exp(θ.V2_idade * (idade - ref.idade))
    V3 = θ.V3 * mm * opi(θ.V3_opioide)
    # depurações (PAS L324–L325, L327, L331)
    CL = (masculino ? θ.CL_m : θ.CL_f) * (peso / ref.peso)^0.75 * mat_cl * opi(θ.CL_opioide)
    Q2 = θ.Q2 * (V2 / θ.V2)^0.75 * (1 + θ.Q2_mat * (1 - q3_imaturo))
    Q3 = θ.Q3 * (V3 / θ.V3)^0.75 * mat_q3
    return (V1 = V1, V2 = V2, V3 = V3, CL = CL, Q2 = Q2, Q3 = Q3)
end

"""
    matrizes_pk(pk; dt = 1.0)

Matrizes da solução exata do modelo de 3 compartimentos (sangue, músculo,
gordura) para um passo de `dt` segundos com a taxa constante no passo, como o
PAS (pk_models.py L630–L664): x(t+dt) = Φ·x(t) + Γ·taxa. A Cp é x[1].
"""
function matrizes_pk(pk; dt = 1.0)
    # taxas de transferência entre compartimentos, de 1/min para 1/s (PAS L526–L536)
    k10, k12, k13 = pk.CL / pk.V1, pk.Q2 / pk.V1, pk.Q3 / pk.V1
    k21, k31 = pk.Q2 / pk.V2, pk.Q3 / pk.V3
    A = [-(k10 + k12 + k13) k12 k13;
         k21 -k21 0.0;
         k31 0.0 -k31] ./ 60
    B = [1 / pk.V1, 0.0, 0.0]          # a dose entra no sangue (compartimento 1)
    Φ = exp(A .* dt)
    Γ = A \ ((Φ - I) * B)
    return Φ, Γ
end

"""
    simular_pk(pk, taxa; dt = 1.0)

Concentração plasmática (µg/ml) a cada `dt` segundos, para uma infusão `taxa`
(mg/s) dada a cada `dt` segundos. `pk` vem de `parametros_eleveld`.
"""
function simular_pk(pk, taxa; dt = 1.0)
    Φ, Γ = matrizes_pk(pk; dt)
    x = zeros(3)
    cp = zeros(length(taxa))
    for k in eachindex(taxa)
        cp[k] = x[1]                   # Cp no início do segundo k
        x = Φ * x .+ Γ .* taxa[k]
    end
    return cp
end

"""
    taxa_bolus(dose_mg, duracao_s, n; dt = 1.0)

Taxa de infusão (mg/s) de um bolus de `dose_mg` injetado em `duracao_s` segundos,
em `n` amostras de `dt` segundos.
"""
taxa_bolus(dose_mg, duracao_s, n; dt = 1.0) = [t < duracao_s ? dose_mg / duracao_s : 0.0 for t in dt .* (0:n-1)]
