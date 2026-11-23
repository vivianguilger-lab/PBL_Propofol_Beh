# =============================================================================
# windkessel.jl — a árvore arterial: Windkessel de 4 elementos
#
# É o centro do modelo: transforma o fluxo que sai do coração, batimento a
# batimento, em pressão arterial (PAS, PAD e PAM). Analogia elétrica (aula de
# Modelagem do Sistema Cardiovascular): pressão = tensão, fluxo = corrente.
#
#            Q(t) →──┬── Rc ──┬──────┬──────┐
#                    └── L  ──┘      │      │
#                                    C      Rp        P(t) = pressão aórtica
#                                    │      │
#            ────────────────────────┴──────┘
#
#   Rp  resistência periférica (arteríolas): define a pressão média
#   C   complacência arterial: guarda sangue na sístole e devolve na diástole
#   Rc  impedância característica da aorta proximal
#   L   inertância do sangue (em paralelo com Rc)
#
# Modelo: Stergiopulos N, Westerhof BE, Westerhof N. Total arterial inertance as
# the fourth element of the windkessel model. Am J Physiol 1999;276:H81–8.
# Implementação de referência: windkessel-id (Pigot et al., licença MIT),
# github.com/hpigot/windkessel-id, já portada pelo grupo no PBL_6
# (src/windkessel.jl); as equações são as mesmas.
#
# Unidades: pressão em mmHg, fluxo em mL/s, volume em mL, tempo em s;
# Rp e Rc em mmHg·s/mL, C em mL/mmHg, L em mmHg·s²/mL.
# =============================================================================

# Parâmetros humanos de Stergiopulos (1999), batimento tipo A (windkessel-id;
# PBL_6 src/windkessel.jl, params_stergiopulos1999total). Rp e C são trocados pelos
# do paciente (parametrizar_windkessel); Rc e L ficam estes.
const WK_STERGIOPULOS = (Rp = 0.79, C = 1.22, Rc = 0.056, L = 0.0051)

# Passo de tempo dentro do batimento: 5 ms, a amostragem do batimento de Stergiopulos.
const H_WK = 0.005

# Forma de onda do fluxo aórtico: o batimento tipo A de Stergiopulos (1999),
# digitalizado pelo windkessel-id (166 amostras a cada 5 ms, cerca de 0,83 s).
# Guardada com média 1, para ser esticada para 60/FC e escalada pelo volume sistólico.
const ARQUIVO_FLUXO = joinpath(@__DIR__, "..", "dados", "stergiopulos1999", "aflow.csv")
const FORMA_FLUXO = let
    q = [parse(Float64, split(l, ",")[2]) for l in eachline(ARQUIVO_FLUXO) if !isempty(strip(l))]
    q ./ (sum(q) / length(q))
end

"""
    fluxo_no_batimento(forma, k, n)

Fluxo relativo (média 1 no batimento) na amostra `k` (0 a n−1) de um batimento com
`n` amostras: a forma de Stergiopulos esticada para o tamanho do batimento, com
interpolação linear (PBL_6 src/acoplamento.jl, fluxo_batimento).
"""
@inline function fluxo_no_batimento(forma, k, n)
    m = length(forma)
    s = k * m / n
    i = floor(Int, s)
    f = s - i
    return (1 - f) * forma[i + 1] + f * forma[mod(i + 1, m) + 1]
end

"""
    matrizes_wk(Rp, C, Rc, L; h = H_WK)

Windkessel de 4 elementos (forma paralela) em espaço de estados, já discretizado
para o passo `h` com a taxa constante no passo (solução exata):

    dx1/dt = −x1/(Rp·C) + Q            x1 = volume guardado em C
    dx2/dt = −(Rc/L)·x2 + Rc·Q         x2 = fluxo no ramo de L, × L
    P      = x1/C − (Rc/L)·x2 + Rc·Q

Como a matriz é diagonal, a discretização é elemento a elemento:
a = e^(λh) e b = (a − 1)/λ · (entrada). Devolve `(a, b, c, d)`.
"""
function matrizes_wk(Rp, C, Rc, L; h = H_WK)
    λ = (-1 / (Rp * C), -Rc / L)
    a = exp.(λ .* h)
    b = (a .- 1) ./ λ .* (1.0, Rc)
    c = (1 / C, -Rc / L)
    return a, b, c, Rc
end

"""
    batimento(Rp, C, Rc, L, fc, vs, x; h = H_WK, forma = FORMA_FLUXO, guardar = false)

Simula UM batimento: o coração ejeta `vs` mL com a forma de Stergiopulos, num
batimento de 60/`fc` s, a partir do estado `x` do Windkessel. Devolve
`(pas, pad, pam, x, duracao)`: pressão máxima, mínima e média do batimento (mmHg),
o estado no fim e a duração (s). Com `guardar = true`, devolve também a pressão e o
fluxo amostra a amostra (`pressao`, `fluxo`).
"""
function batimento(Rp, C, Rc, L, fc, vs, x; h = H_WK, forma = FORMA_FLUXO, guardar = false)
    n = max(round(Int, 60 / fc / h), 2)
    # mL/s por unidade da forma, para que o fluxo médio do batimento seja exatamente o débito
    # VS × FC/60: a duração do batimento é arredondada para múltiplos de h (n·h ≈ 60/FC)
    escala = vs * fc / 60 * n / sum(fluxo_no_batimento(forma, k, n) for k in 0:n-1)
    a, b, c, d = matrizes_wk(Rp, C, Rc, L; h)
    x1, x2 = x
    pmax, pmin, soma = -Inf, Inf, 0.0
    pressao = guardar ? zeros(n) : Float64[]
    fluxo = guardar ? zeros(n) : Float64[]
    for k in 0:n-1
        q = escala * fluxo_no_batimento(forma, k, n)
        p = c[1] * x1 + c[2] * x2 + d * q
        pmax, pmin, soma = max(pmax, p), min(pmin, p), soma + p
        if guardar
            pressao[k + 1], fluxo[k + 1] = p, q
        end
        x1 = a[1] * x1 + b[1] * q
        x2 = a[2] * x2 + b[2] * q
    end
    r = (pas = pmax, pad = pmin, pam = soma / n, x = (x1, x2), duracao = n * h)
    return guardar ? merge(r, (pressao = pressao, fluxo = fluxo)) : r
end

"""
    estado_periodico(Rp, C, Rc, L, fc, vs; h = H_WK)

Estado inicial do batimento em regime permanente (o batimento se repete igual):
partindo do repouso, x_n = a^n·x0 + x_n⁰, e x0 = x_n pede x0 = x_n⁰/(1 − a^n)
(windkessel-id, solvex0; PBL_6 src/windkessel.jl).
"""
function estado_periodico(Rp, C, Rc, L, fc, vs; h = H_WK)
    n = max(round(Int, 60 / fc / h), 2)
    a, = matrizes_wk(Rp, C, Rc, L; h)
    xn = batimento(Rp, C, Rc, L, fc, vs, (0.0, 0.0); h).x
    return (xn[1] / (1 - a[1]^n), xn[2] / (1 - a[2]^n))
end

"""
    parametrizar_windkessel(; pam, pad, fc, dc = nothing, peso, altura,
                            Rc = WK_STERGIOPULOS.Rc, L = WK_STERGIOPULOS.L)

O Windkessel do paciente (o "gêmeo digital" da árvore arterial), a partir do monitor:

- débito cardíaco basal DC0: o medido (L/min) ou, se não houver, 3,0 L/min/m² ×
  superfície corporal de Mosteller (como no PBL_6: índice cardíaco normal em repouso);
- volume sistólico basal VS0 = DC0/FC0;
- Rp = PAM/DC0: no regime médio, C não deixa passar fluxo e L curto-circuita Rc,
  então PAM = Rp × fluxo médio (lei de Ohm da circulação);
- C: identificada para que a PAD do batimento simulado seja a medida, por bissecção
  (com a PAM fixa, a PAD sobe quando C aumenta). Usa-se a PAD, e não a PAS, porque o
  Windkessel dá a pressão central (aórtica): a PAM e a PAD quase não mudam da aorta
  ao braço, mas a PAS do manguito é amplificada na periferia (Agabiti-Rosei et al.,
  2007), então a PAS central simulada fica abaixo da PAS medida no braço;
- Rc e L: de Stergiopulos (1999).

Devolve `(Rp, C, Rc, L, dc, vs, pas, pad, pam, x0)`, com a PAS, a PAD e a PAM do
batimento basal simulado (centrais) e o estado periódico inicial.
"""
function parametrizar_windkessel(; pam, pad, fc, dc = nothing, peso, altura,
                                 Rc = WK_STERGIOPULOS.Rc, L = WK_STERGIOPULOS.L)
    pad < pam || error("A PAD deve ser menor que a PAM.")
    dc0 = dc === nothing ? 3.0 * sqrt(altura * peso / 3600) : dc     # L/min
    q0 = dc0 * 1000 / 60                                               # mL/s
    vs0 = q0 * 60 / fc                                                 # mL
    Rp = pam / q0
    basal(C) = batimento(Rp, C, Rc, L, fc, vs0, estado_periodico(Rp, C, Rc, L, fc, vs0))
    lo, hi = 0.05, 20.0                                                # mL/mmHg
    if basal(hi).pad < pad         # nem com C muito grande a PAD sobe tanto: fica no limite
        C = hi
    elseif basal(lo).pad > pad     # nem com C muito pequena a PAD desce tanto
        C = lo
    else
        for _ in 1:60
            m = sqrt(lo * hi)
            basal(m).pad < pad ? (lo = m) : (hi = m)
        end
        C = sqrt(lo * hi)
    end
    x0 = estado_periodico(Rp, C, Rc, L, fc, vs0)
    b = batimento(Rp, C, Rc, L, fc, vs0, x0)
    return (Rp = Rp, C = C, Rc = Rc, L = L, dc = dc0, vs = vs0, pas = b.pas, pad = b.pad, pam = b.pam, x0 = x0)
end
