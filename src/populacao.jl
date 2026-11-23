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
#
# ATUALIZAÇÃO (remifentanil): também variam os 6 parâmetros da farmacocinética de Minto
# e os 3 parâmetros do remifentanil no Su que têm variabilidade publicada (ver abaixo).
# =============================================================================

# Variâncias ω² entre pessoas (log-normais)
const ω²_ELEVELD = (V1 = 0.610, V2 = 0.565, V3 = 0.597, CL = 0.265, Q2 = 0.346, Q3 = 0.209)  # PAS L346–L351
const ω²_SU_C50_RPT = 0.44                                                                     # PAS L1009
# Variabilidade do BIS (Eleveld 2018): Ce50 (o PAS guarda como coeficiente de variação 0,523 e converte
# em ω² = ln(1 + 0,523²) = 0,242; pd_models.py L223 e L285) e ke0 (pk_models.py L352).
const ω²_ELEVELD_BIS = (ce50 = log(1 + 0.523^2), ke0 = 0.702)
# Variabilidade do remifentanil no Su 2023. No PAS, estas três são ADITIVAS (valor + η, e não × exp(η)):
#   efeito máximo do remi na RPT: η ~ Normal(0, variância 0,449)            pd_models.py L1010, L1077
#   inclinações do remi na FC e no VS: η sorteados juntos (normal bivariada) pd_models.py L1007, L1079–L1081
const ω²_SU_EMAX_REMI_RPT = 0.449
const Ω_SU_REMI_FC_VS = [0.00382 0.00329; 0.00329 0.00868]    # covariância (FC, VS)

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
    sortear_remifentanil(pr, psu, rng; ω²_pr = ω²_MINTO)

Parte do remifentanil de uma versão do paciente: a farmacocinética de Minto `pr`
(× exp(η), como a do propofol) e os 3 parâmetros do remifentanil de `psu` com
variabilidade publicada (somados a η, como no PAS).
"""
function sortear_remifentanil(pr, psu, rng; ω²_pr = ω²_MINTO)
    η(ω²) = sqrt(ω²) * randn(rng)
    pr_versao = (V1 = pr.V1 * exp(η(ω²_pr.V1)), V2 = pr.V2 * exp(η(ω²_pr.V2)),
                 V3 = pr.V3 * exp(η(ω²_pr.V3)), CL = pr.CL * exp(η(ω²_pr.CL)),
                 Q2 = pr.Q2 * exp(η(ω²_pr.Q2)), Q3 = pr.Q3 * exp(η(ω²_pr.Q3)))
    # FC e VS sorteados juntos: η = L·z, com L·Lᵀ = covariância (fatoração de Cholesky)
    η_fc, η_vs = cholesky(Ω_SU_REMI_FC_VS).L * randn(rng, 2)
    psu_versao = trocar(trocar(trocar(psu, :emax_remi_rpt, psu.emax_remi_rpt + η(ω²_SU_EMAX_REMI_RPT)),
                               :sl_remi_fc, psu.sl_remi_fc + η_fc), :sl_remi_vs, psu.sl_remi_vs + η_vs)
    return pr_versao, psu_versao
end

"""
    simular_populacao(idade, peso, altura, sexo, pam_basal, fc_basal, taxa; wk,
                      taxa_remi = nothing, antes_s = 0, n = 500, semente = 1,
                      ω²_pk = ω²_ELEVELD, ω²_c50 = ω²_SU_C50_RPT, ω²_bis = ω²_ELEVELD_BIS)

Simula `n` versões do paciente com a mesma infusão de propofol `taxa` (mg/s, uma
por segundo, de 0 em diante) e, se houver, de remifentanil `taxa_remi` (µg/s, uma
por segundo, começando `antes_s` segundos antes do propofol). Devolve o tempo `t`
(s, 0 = início do propofol) e, uma curva por versão: `curvas` (PAM, mmHg), `pas` e
`pad` (centrais, mmHg), `fc` (bpm), `dc` (débito cardíaco, L/min), `cp` (Cp do
propofol, µg/ml), `cr` (Cp do remifentanil, ng/ml), `rpt` e `vs` (resistência e
volume sistólico, como fração do basal) e `bis` (BIS no monitor).

A pressão de cada versão sai do Windkessel do paciente `wk` (parametrizar_windkessel),
em malha fechada com o modelo de Su (simular_circulacao). A `semente` fixa o sorteio:
a mesma semente dá sempre os mesmos pacientes, com ou sem remifentanil.
"""
function simular_populacao(idade, peso, altura, sexo, pam_basal, fc_basal, taxa; wk,
                           taxa_remi = nothing, antes_s = 0, n = 500, semente = 1,
                           ω²_pk = ω²_ELEVELD, ω²_c50 = ω²_SU_C50_RPT, ω²_bis = ω²_ELEVELD_BIS)
    rng = Xoshiro(semente)
    com_remi = taxa_remi !== nothing && any(>(0), taxa_remi)
    antes = com_remi ? antes_s : 0
    com_remi && length(taxa_remi) != antes + length(taxa) &&
        error("taxa_remi deve ir de -antes_s até o fim da simulação.")
    # com remifentanil, a farmacocinética do propofol usa os termos de Eleveld para opioide junto
    pk = parametros_eleveld(idade, peso, altura, sexo; opioide = com_remi)
    pr = parametros_minto(idade, peso, altura, sexo)
    psu = ParametrosSu()
    pb = ParametrosBIS()
    tp = vcat(zeros(antes), taxa)          # propofol: zero enquanto só o remifentanil corre
    fica = antes + 1:antes + length(taxa)  # o que é mostrado: do início do propofol em diante
    versoes = map(1:n) do _
        pk_versao, psu_versao = sortear_paciente(pk, psu, rng; ω²_pk, ω²_c50)
        # BIS da versão: Ce50 e ke0 também sorteados (sorteados depois dos outros, para não mudar a PAM)
        pb_versao = trocar(trocar(pb, :ce50_35, pb.ce50_35 * exp(sqrt(ω²_bis.ce50) * randn(rng))),
                           :ke0_70, pb.ke0_70 * exp(sqrt(ω²_bis.ke0) * randn(rng)))
        # remifentanil da versão: sorteado sempre, para que cada versão seja a mesma com ou sem remifentanil
        pr_versao, psu_versao = sortear_remifentanil(pr, psu_versao, rng)
        cp = simular_pk(pk_versao, tp)
        cr = com_remi ? simular_pk(pr_versao, taxa_remi) : zeros(length(tp))
        # circulação: Windkessel do paciente, com Rp, VS e FC dados pelo modelo de Su desta versão
        h = simular_circulacao(cp, cr, idade, wk; pam0 = pam_basal, fc0 = fc_basal, p = psu_versao)
        (cp = cp[fica], cr = cr[fica], pam = h.pam[fica], pas = h.pas[fica], pad = h.pad[fica], dc = h.dc[fica],
         fc = h.fc[fica], rpt = h.rpt_rel[fica], vs = h.vs_rel[fica],
         bis = simular_bis(pk_versao, tp, idade, peso; p = pb_versao).bis[fica])
    end
    return (t = collect(0.0:length(taxa) - 1),
            curvas = [v.pam for v in versoes],        # PAM (mmHg)
            pas = [v.pas for v in versoes],           # PAS central (mmHg)
            pad = [v.pad for v in versoes],           # PAD (mmHg)
            dc = [v.dc for v in versoes],             # débito cardíaco (L/min)
            fc = [v.fc for v in versoes],             # FC (bpm)
            cp = [v.cp for v in versoes],             # Cp do propofol (µg/ml)
            cr = [v.cr for v in versoes],             # Cp do remifentanil (ng/ml)
            rpt = [v.rpt for v in versoes],           # RPT (fração do basal)
            vs = [v.vs for v in versoes],             # VS (fração do basal)
            bis = [v.bis for v in versoes])           # BIS (0–100)
end
