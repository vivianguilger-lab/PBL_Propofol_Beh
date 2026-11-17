# =============================================================================
# populacao.jl — pacientes virtuais: várias "versões" do mesmo paciente
#
# Cada versão tem os parâmetros sorteados dentro da variabilidade entre pessoas
# publicada nos artigos. A forma é a dos próprios artigos (log-normal):
#     valor da versão = valor típico × exp(η),   η ~ Normal(média 0, variância ω²)
#
# Variabilidade usada (ω²):
#   Eleveld 2018: V1, V2, V3, CL, Q2 e Q3 (PAS pk_models.py L346–L351).
#   Su 2023: C50 do efeito na RPT (PAS pd_models.py L1009). É a única variabilidade
#   do efeito do propofol que o PAS traz; os outros parâmetros do Su ficam fixos.
#   A variabilidade dos valores basais (RPT, VS e FC) do Su não entra, porque
#   aqui os basais são os medidos no paciente.
#
# O sorteio segue a definição acima, e não o código do PAS para o C50 da RPT:
# em pd_models.py L1054 e L1075, o PAS passa ω como média do sorteio, e não como
# desvio-padrão (na farmacocinética, L518, ele faz certo).7
#
# PONTO DE ATENÇÃO, NÃO SEI SE TALVEZ SEJA MELHOR MUDAR ISSO:
# Varia: os 6 parâmetros da farmacocinética (Eleveld) e o C50 da resistência (Su). É toda a variabilidade publicada que o PAS traz para o efeito do propofol.
# Não varia: os outros 11 parâmetros do Su, porque não há variabilidade publicada no PAS para eles.
# Consequência: a faixa provavelmente é mais estreita do que a realidade. É exatamente o que o VitalDB vai corrigir depois.
# =============================================================================

# Variâncias ω² entre pessoas (log-normais)
const ω²_ELEVELD = (V1 = 0.610, V2 = 0.565, V3 = 0.597, CL = 0.265, Q2 = 0.346, Q3 = 0.209)  # PAS L346–L351
const ω²_SU_C50_RPT = 0.44                                                                     # PAS L1009

# Cópia de um registro (como ParametrosSu) com um campo trocado.
trocar(p, campo, valor) = typeof(p)((f == campo ? valor : getfield(p, f) for f in fieldnames(typeof(p)))...)

"""
    sortear_paciente(pk, psu, rng; ω²_pk = ω²_ELEVELD, ω²_c50 = ω²_SU_C50_RPT)

Uma versão do paciente: os parâmetros farmacocinéticos `pk` (de
`parametros_eleveld`) e o C50 da RPT de `psu` (`ParametrosSu`), cada um
multiplicado por exp(η) sorteado. `rng` é o gerador de números aleatórios.
"""
function sortear_paciente(pk, psu, rng; ω²_pk = ω²_ELEVELD, ω²_c50 = ω²_SU_C50_RPT)
    # η: número sorteado de uma normal com média 0 e variância ω²
    η(ω²) = sqrt(ω²) * randn(rng)
    pk_versao = (V1 = pk.V1 * exp(η(ω²_pk.V1)), V2 = pk.V2 * exp(η(ω²_pk.V2)),
                 V3 = pk.V3 * exp(η(ω²_pk.V3)), CL = pk.CL * exp(η(ω²_pk.CL)),
                 Q2 = pk.Q2 * exp(η(ω²_pk.Q2)), Q3 = pk.Q3 * exp(η(ω²_pk.Q3)))
    psu_versao = trocar(psu, :c50_rpt, psu.c50_rpt * exp(η(ω²_c50)))
    return pk_versao, psu_versao
end

"""
    simular_populacao(idade, peso, altura, sexo, pam_basal, fc_basal, taxa;
                      n = 500, semente = 1, ω²_pk = ω²_ELEVELD, ω²_c50 = ω²_SU_C50_RPT)

Simula `n` versões do paciente com a mesma infusão `taxa` (mg/s, uma por
segundo). Devolve `(t, curvas)`: o tempo (s) e uma curva de PAM (mmHg) por
versão. A `semente` fixa o sorteio: a mesma semente dá sempre as mesmas curvas.
"""
function simular_populacao(idade, peso, altura, sexo, pam_basal, fc_basal, taxa;
                           n = 500, semente = 1, ω²_pk = ω²_ELEVELD, ω²_c50 = ω²_SU_C50_RPT)
    rng = Xoshiro(semente)
    pk = parametros_eleveld(idade, peso, altura, sexo)
    psu = ParametrosSu()
    curvas = map(1:n) do _
        pk_versao, psu_versao = sortear_paciente(pk, psu, rng; ω²_pk, ω²_c50)
        cp = simular_pk(pk_versao, taxa)
        simular_hemodinamica(cp, pam_basal, fc_basal, idade; p = psu_versao).pam
    end
    return (t = collect(0.0:length(taxa) - 1), curvas = curvas)
end
