# =============================================================================
# alvo.jl — entrada por concentração-alvo (como uma bomba de TCI)
#
# O médico escolhe a concentração plasmática de propofol que quer atingir
# (µg/ml) e esta parte calcula, segundo a segundo, a taxa de infusão (mg/s) que
# leva o modelo farmacocinético exatamente a essa concentração e a mantém.
#
# Escolhas:
#   - alvo no plasma (Cp), e não no sítio efetor: o modelo de Su 2023 usa a Cp,
#     e o alvo no sítio efetor exigiria um ke0 (docs/decisoes.md);
#   - a taxa é calculada com o paciente TÍPICO, como faz uma bomba real; cada
#     versão do paciente (populacao.jl) recebe essa mesma taxa e, por ter outra
#     farmacocinética, chega a uma Cp diferente do alvo;
#   - sem limite de velocidade da bomba por padrão (`taxa_max = Inf`): o limite
#     depende do equipamento e pode ser informado.
# =============================================================================

"""
    taxa_alvo(pk, alvo, n; dt = 1.0, taxa_max = Inf)

Taxa de infusão (mg/s), uma a cada `dt` segundos durante `n` amostras, para que
a Cp do modelo `pk` chegue ao `alvo` (µg/ml) e fique nele. `alvo` pode ser um
número (fixo) ou uma lista com um alvo por amostra (para mudar o alvo no tempo).
A taxa nunca é negativa: se a Cp está acima do alvo, a bomba para até ela cair.
"""
function taxa_alvo(pk, alvo, n; dt = 1.0, taxa_max = Inf)
    Φ, Γ = matrizes_pk(pk; dt)
    alvos = alvo isa Number ? fill(float(alvo), n) : alvo
    x = zeros(3)
    taxa = zeros(n)
    for k in 1:n
        # Cp que haveria no próximo passo sem infusão: (Φ·x)[1]. Cada mg/s a mais
        # soma Γ[1] µg/ml. A taxa certa é a que completa o que falta para o alvo.
        sem_infusao = (Φ * x)[1]
        taxa[k] = clamp((alvos[k] - sem_infusao) / Γ[1], 0.0, taxa_max)
        x = Φ * x .+ Γ .* taxa[k]
    end
    return taxa
end
