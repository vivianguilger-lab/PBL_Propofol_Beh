# =============================================================================
# circulacao.jl — a circulação em malha fechada: coração + Windkessel + barorreflexo
#
#   coração (fonte de fluxo) ── Q(t) ──▶ Windkessel de 4 elementos ──▶ PAS, PAD, PAM
#        ▲  VS e FC                          ▲  Rp                         │
#        │                                   │                             │
#        └───── Su 2023: fármacos + barorreflexo (feedback da PAM) ◀───────┘
#
# A cada batimento:
#   1. o estado do modelo de Su dá Rp = Rp0 × RPT, VS = VS0 × VS_rel e FC = FC0 × FC_rel
#      (frações do basal; Rp0, VS0 e FC0 são os do paciente, windkessel.jl);
#   2. o Windkessel simula o batimento (passo de 5 ms) e dá PAS, PAD e PAM;
#   3. a PAM média do batimento volta ao modelo de Su como feedback do barorreflexo,
#      e o estado de Su avança a duração do batimento (Runge-Kutta de 4ª ordem),
#      com as concentrações de propofol e remifentanil daquele momento.
#
# Congruência: em regime, a PAM do Windkessel é Rp × fluxo médio = (PAM0/DC0)·RPT ×
# DC0·VS_rel·FC_rel = PAM0·RPT·VS_rel·FC_rel, exatamente a PAM do modelo de Su. O
# Windkessel acrescenta a onda de pressão de cada batimento (PAS e PAD) e a dinâmica
# da complacência; a PAM só difere da do Su nos transitórios (teste em runtests.jl).
# =============================================================================

"""
    simular_circulacao(cp, cr, idade, wk; pam0, fc0, p = ParametrosSu())

Simula a circulação batimento a batimento.

- `cp`, `cr`: Cp do propofol (µg/ml) e do remifentanil (ng/ml), uma por segundo,
  do segundo 0 ao fim (`cr = nothing`: sem remifentanil).
- `wk`: o Windkessel do paciente (`parametrizar_windkessel`); `pam0`, `fc0`: PAM e FC basais.
- `p`: parâmetros do modelo de Su (podem ser os de um paciente virtual).

Devolve, uma por segundo: `pam`, `pas`, `pad` (mmHg, do batimento em curso), `fc`
(bpm), `vs_rel`, `rpt_rel` (frações do basal) e `dc` (débito cardíaco, L/min).
"""
function simular_circulacao(cp, cr, idade, wk; pam0, fc0, p = ParametrosSu())
    n = length(cp)
    r = cr === nothing ? zeros(n) : cr
    conc(v, t) = (i = clamp(floor(Int, t), 0, n - 2); f = clamp(t - i, 0.0, 1.0); (1 - f) * v[i + 1] + f * v[i + 2])
    y = estado_basal(p)
    x = wk.x0
    pam_rel = 1.0
    saida = (pam = zeros(n), pas = zeros(n), pad = zeros(n), fc = zeros(n), vs_rel = zeros(n), rpt_rel = zeros(n), dc = zeros(n))
    t, s = 0.0, 0                     # tempo (s) e próximo segundo a preencher
    while s < n
        a = aparentes(y, p)
        fc, vs, Rp = fc0 * a.fc, wk.vs * a.vs, wk.Rp * a.rpt
        b = batimento(Rp, wk.C, wk.Rc, wk.L, fc, vs, x)
        # os segundos que caem neste batimento recebem os valores dele
        while s < n && s < t + b.duracao
            saida.pam[s + 1], saida.pas[s + 1], saida.pad[s + 1] = b.pam, b.pas, b.pad
            saida.fc[s + 1], saida.vs_rel[s + 1], saida.rpt_rel[s + 1] = fc, a.vs, a.rpt
            saida.dc[s + 1] = vs * fc / 1000
            s += 1
        end
        # barorreflexo: a PAM deste batimento volta ao modelo de Su
        pam_rel = b.pam / pam0
        tm = t + b.duracao / 2
        y = passo_rk4(y, conc(cp, tm), conc(r, tm), idade, b.duracao, p; pam_rel)
        x = b.x
        t += b.duracao
    end
    return saida
end

"""
    onda_de_pressao(cp, cr, idade, wk; pam0, fc0, momentos, p = ParametrosSu())

Para o gráfico da onda: simula como `simular_circulacao` e guarda, amostra a amostra
(5 ms), a pressão e o fluxo do batimento em curso em cada instante de `momentos` (s).
Devolve uma lista de `(t, pressao, fluxo, pas, pad, pam, fc)`.
"""
function onda_de_pressao(cp, cr, idade, wk; pam0, fc0, momentos, p = ParametrosSu())
    n = length(cp)
    r = cr === nothing ? zeros(n) : cr
    conc(v, t) = (i = clamp(floor(Int, t), 0, n - 2); f = clamp(t - i, 0.0, 1.0); (1 - f) * v[i + 1] + f * v[i + 2])
    y, x, t = estado_basal(p), wk.x0, 0.0
    faltam = sort(collect(momentos))
    ondas = []
    while !isempty(faltam) && t < n
        a = aparentes(y, p)
        fc = fc0 * a.fc
        b = batimento(wk.Rp * a.rpt, wk.C, wk.Rc, wk.L, fc, wk.vs * a.vs, x; guardar = true)
        while !isempty(faltam) && faltam[1] < t + b.duracao
            push!(ondas, (t = popfirst!(faltam), pressao = b.pressao, fluxo = b.fluxo, pas = b.pas, pad = b.pad, pam = b.pam, fc = fc))
        end
        tm = t + b.duracao / 2
        y = passo_rk4(y, conc(cp, tm), conc(r, tm), idade, b.duracao, p; pam_rel = b.pam / pam0)
        x, t = b.x, t + b.duracao
    end
    return ondas
end
