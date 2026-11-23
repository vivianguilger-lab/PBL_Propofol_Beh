# =============================================================================
# hemodinamica.jl — efeito do propofol sobre a circulação
#
# Modelo: Su H, Koomen JV, Eleveld DJ, Struys MMRF, Colin PJ. Pharmacodynamic
# mechanism-based interaction model for the haemodynamic effects of remifentanil
# and propofol in healthy volunteers. Br J Anaesth 2023;131(2):222–33.
# doi:10.1016/j.bja.2023.04.043
#
# Só a parte do propofol: com remifentanil = 0, todos os termos do remifentanil                  DÚVIDA: Eu não quero medir o remifetanil, 
# do modelo original valem exatamente zero.                                                      pq nosso foco é o propofol, mas se eu igualo 
#                                                                                                ele à 0 não pode mudar a análise.
# Implementação de referência: Python Anesthesia Simulator (PAS), commit 4afc03a,                Verificação: Zerar todos os coeficientes do 
# src/python_anesthesia_simulator/pd_models.py, classe Hemo_meca_PD_model (L890).                modelo invalida alguma hipotese dele
#                                                                                                Ou seja, se zero pode fazer parte do conjunto 
#                                                                                                solução dos coeficientes do remifetanil, de 
#                                                                                                acordo com modelo escolhido
#
#
#
#                                                                                                
#                                                                                                DÚVIDA: acho q os 21 casos do VitalDB selecionados no PBL_6
#                                                                                                estão todos com coadministração ativa (coadmin_status = active), 
#                                                                                                ou seja, com remifentanil junto. Vamos usar esses casos? 
#                                                                                                Se sim, acho q o remifentanil precisa entrar no modelo
#
# ATUALIZAÇÃO (remifentanil): o remifentanil agora entra no modelo, com os coeficientes
# publicados por Su 2023 (PAS L971–L978) e a Cp do remifentanil de Minto (remifentanil.jl).
# Com Cp do remifentanil = 0, tudo fica exatamente igual à versão só com propofol.
#
# ATUALIZAÇÃO (Windkessel): o modelo central é o Windkessel de 4 elementos (windkessel.jl).
# Este arquivo dá a ele os parâmetros que mudam com os fármacos: Rp (= RPT), VS e FC, com
# os valores publicados por Su 2023 no lugar dos valores ilustrativos do modelo anterior
# (Imax, Ce50, γ e o fator de sensibilidade s). Na simulação acoplada (circulacao.jl), o
# feedback do barorreflexo usa a PAM calculada pelo Windkessel (argumento pam_rel).
# =============================================================================
# Curva sigmoide (Hill): vai de 0 (sem propofol) a 1 (efeito máximo).                           
# c50 = concentração que dá metade do efeito; γ = inclinação da curva.
# PAS: pd_models.py L11.
fsig(x, c50, γ) = x^γ / (c50^γ + x^γ)

# Parâmetros publicados por Su 2023, como estão no PAS (pd_models.py L958–L969).
# As taxas estão em 1/min no artigo; aqui são divididas por 60 (1/s).
Base.@kwdef struct ParametrosSu
    k_out::Float64      = 0.072 / 60  # renovação de RPT, VS e FC (1/s)            L958              As linhas de onde esses valores foram retirados do
    fb::Float64         = -0.661      # força do feedback da PAM (expoente)        L959              repositório estão referenciados, mas estou confusa se
    α::Float64          = 0.312       # quanto a FC altera o VS                    L960              é melhor deixar assim ou substituir por valores baseados
    k_ltde::Float64     = 0.067 / 60  # velocidade do efeito dependente do tempo   L961              no banco de dados
    ltde_vs::Float64    = 0.0899      # efeito dependente do tempo no VS (fração)  L962
    ltde_fc::Float64    = 0.121       # efeito dependente do tempo na FC (fração)  L963
    c50_rpt::Float64    = 3.21        # Cp com metade do efeito na RPT (µg/ml)     L964
    emax_rpt::Float64   = -0.778      # efeito máximo na RPT (−77,8%)              L965
    γ_rpt::Float64      = 1.83        # inclinação do efeito na RPT                L966
    c50_vs::Float64     = 0.44        # Cp com metade do efeito no VS (µg/ml)      L967
    emax_vs35::Float64  = -0.154      # efeito máximo no VS aos 35 anos (−15,4%)   L968
    k_idade_vs::Float64 = 0.0333      # aumento do efeito no VS por ano de idade   L969
    # remifentanil (Cp em ng/ml) e interação com o propofol
    c50_remi_rpt::Float64  = 4.59     # Cp do remi com metade do efeito na RPT (ng/ml)   L971
    emax_remi_rpt::Float64 = -1.0     # efeito máximo do remi na RPT                   L972
    sl_remi_fc::Float64    = 0.0327   # efeito do remi na FC, por ng/ml                L973
    sl_remi_vs::Float64    = 0.0581   # efeito do remi no VS, por ng/ml                L974
    int_fc::Float64        = -0.119   # interação na FC: o propofol muda o efeito do remi      L975
    c50_int_fc::Float64    = 0.20     # Cp do propofol com metade dessa interação (µg/ml)    L976
    int_rpt::Float64       = 1.0      # interação na RPT: o remi muda o efeito do propofol     L977
    int_vs::Float64        = -0.212   # interação no VS: o propofol muda o efeito do remi      L978
end

# Efeito máximo do propofol no VS para a idade do paciente (PAS L969):
# −0,154 aos 35 anos, crescendo 3,33% (multiplicativo) por ano acima disso.
emax_vs(p::ParametrosSu, idade) = p.emax_vs35 * exp(p.k_idade_vs * (idade - 35))

# Estado do modelo, sempre como fração do valor basal (1 = basal):
#   y[1] = RPT  (resistência periférica total)
#   y[2] = VS   (volume sistólico, sem o termo dependente do tempo)
#   y[3] = FC   (frequência cardíaca, sem o termo dependente do tempo)
#   y[4] = termo dependente do tempo do VS
#   y[5] = termo dependente do tempo da FC
# Trabalhar em frações faz o resultado de PAM e FC não depender do VS basal
# (que não é medido na rotina).
estado_basal(p::ParametrosSu) = (1.0, 1.0, 1.0, p.ltde_vs, p.ltde_fc)

# Valores "aparentes" (os que seriam medidos), como fração do basal medido.
# FC aparente = FC + termo do tempo (PAS L1166); VS aparente = (VS + termo do
# tempo) × (1 − α·ln(FC aparente / FC basal aparente)) (PAS L1168);
# PAM = RPT × VS × FC (PAS L1169), dividida pela PAM basal (PAS L1171).
function aparentes(y, p::ParametrosSu)
    fc = (y[3] + y[5]) / (1 + p.ltde_fc)
    vs = (y[2] + y[4]) * (1 - p.α * log(fc)) / (1 + p.ltde_vs)
    rpt = y[1]
    return (pam = rpt * vs * fc, fc = fc, vs = vs, rpt = rpt)
end

# Derivadas do estado (PAS continuous_dynamic, L1118–L1195).
# Cada variável é "produzida" a uma taxa k_out × (PAM relativa)^fb × (1 + efeito
# do propofol) e "removida" a uma taxa k_out × valor atual × (1 − efeito do remi).
# Quando a PAM cai, (PAM relativa)^fb sobe (fb < 0): é o feedback que puxa a PAM de volta.
# cp: Cp do propofol (µg/ml); cr: Cp do remifentanil (ng/ml).
# pam_rel: PAM/PAM basal usada no feedback. Sem ela (nothing), é a do próprio Su (RPT × VS × FC);
# acoplado ao Windkessel, é a PAM média do último batimento calculada pelo Windkessel.
function derivadas(y, cp, cr, idade, p::ParametrosSu; pam_rel = nothing)
    lim(x) = clamp(x, -0.999, 0.999)                                          # L1157–L1161
    # propofol (o remi muda o efeito do propofol na RPT: interação)            L1149–L1150, L1153
    e_rpt = lim((p.emax_rpt + p.int_rpt * fsig(cr, p.c50_remi_rpt, 1)) * fsig(cp, p.c50_rpt, p.γ_rpt))
    e_vs = lim(emax_vs(p, idade) * fsig(cp, p.c50_vs, 1))
    # remifentanil (o propofol muda o efeito do remi no VS e na FC: interação)  L1151, L1152, L1154
    r_rpt = lim(p.emax_remi_rpt * fsig(cr, p.c50_remi_rpt, 1))
    r_vs = lim((p.sl_remi_vs + p.int_vs * fsig(cp, p.c50_vs, 1)) * cr)
    r_fc = lim((p.sl_remi_fc + p.int_fc * fsig(cp, p.c50_int_fc, 1)) * cr)
    retro = (pam_rel === nothing ? aparentes(y, p).pam : pam_rel)^p.fb       # L1171, L1175
    # o termo dependente do tempo só começa a cair depois que o propofol entra (L1184)
    k_t = cp > 0 ? p.k_ltde : 0.0
    return (p.k_out * (retro * (1 + e_rpt) - y[1] * (1 - r_rpt)),   # RPT  L1175–L1176
            p.k_out * (retro * (1 + e_vs) - y[2] * (1 - r_vs)),     # VS   L1179
            p.k_out * (retro - y[3] * (1 - r_fc)),                  # FC   L1182
            -k_t * y[4],                                            # termo do tempo, VS  L1185
            -k_t * y[5])                                            # termo do tempo, FC  L1186
end

# Um passo de Runge-Kutta de 4ª ordem com as Cp constantes no passo.
function passo_rk4(y, cp, cr, idade, dt, p; pam_rel = nothing)
    k1 = derivadas(y, cp, cr, idade, p; pam_rel)
    k2 = derivadas(y .+ dt / 2 .* k1, cp, cr, idade, p; pam_rel)
    k3 = derivadas(y .+ dt / 2 .* k2, cp, cr, idade, p; pam_rel)
    k4 = derivadas(y .+ dt .* k3, cp, cr, idade, p; pam_rel)
    return y .+ dt / 6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end

"""
    simular_hemodinamica(cp, pam_basal, fc_basal, idade; cr = nothing, dt = 1.0, p = ParametrosSu(), subpassos = 10)

Resposta da circulação às concentrações plasmáticas de propofol e remifentanil.

- `cp`: Cp do propofol (µg/ml) a cada `dt` segundos, começando em t = 0.
- `cr`: Cp do remifentanil (ng/ml), nos mesmos instantes (`nothing` = sem remifentanil).
- `pam_basal`, `fc_basal`: PAM (mmHg) e FC (bpm) medidas antes de qualquer fármaco.
- `idade`: anos (muda o efeito do propofol no VS).

Devolve `(t, pam, fc, vs_rel, rpt_rel)`: tempo (s), PAM (mmHg), FC (bpm) e VS e
RPT como fração do basal. Entre duas amostras, usa a média das duas Cp, como o
PAS (pd_models.py L1262).
"""
function simular_hemodinamica(cp, pam_basal, fc_basal, idade; cr = nothing, dt = 1.0, p = ParametrosSu(), subpassos = 10)
    n = length(cp)
    r = cr === nothing ? zeros(n) : cr
    y = estado_basal(p)
    pam, fc, vs, rpt = zeros(n), zeros(n), zeros(n), zeros(n)
    for k in 1:n
        if k > 1
            c = (cp[k - 1] + cp[k]) / 2
            c_r = (r[k - 1] + r[k]) / 2
            for _ in 1:subpassos
                y = passo_rk4(y, c, c_r, idade, dt / subpassos, p)
            end
        end
        a = aparentes(y, p)
        pam[k], fc[k], vs[k], rpt[k] = pam_basal * a.pam, fc_basal * a.fc, a.vs, a.rpt
    end
    return (t = dt .* (0:n-1), pam = pam, fc = fc, vs_rel = vs, rpt_rel = rpt)
end
