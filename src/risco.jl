# =============================================================================
# risco.jl — faixa da PAM e períodos de risco de hipotensão
#
# Recebe várias curvas de PAM possíveis para o mesmo paciente (uma por
# "versão" do paciente) e devolve:
#   1. a faixa da PAM a cada segundo (percentis entre as curvas);
#   2. os períodos de tempo em que a faixa fica em cada nível de risco.
# Não depende de onde vêm as curvas (variabilidade publicada ou VitalDB).
# =============================================================================

# Limiares de PAM (mmHg) que separam os níveis de risco. São definições clínicas:
#   65 — Salmasi V et al. Anesthesiology 2017;126:47–65: PAM abaixo de 65 mmHg
#        associada a lesão renal aguda e lesão miocárdica após cirurgia.
#   55 — Walsh M et al. Anesthesiology 2013;119:507–15: risco maior de lesão renal
#        e miocárdica com PAM abaixo de 55 mmHg.
const LIMIARES = (moderado = 65.0, alto = 55.0)

# Nível de risco de um valor de PAM: baixo (≥ 65), moderado (55 a 65) ou alto (< 55).
nivel(pam; limiares = LIMIARES) =
    pam >= limiares.moderado ? :baixo : pam >= limiares.alto ? :moderado : :alto

# Percentil p (de 0 a 1) de uma lista de números: ordena e interpola em linha
# reta entre os dois valores vizinhos. É o mesmo método do PERCENTIL do Excel
# (PERCENTIL.INC) e do pacote Statistics do Julia.
function percentil(x, p)
    v = sort(x)
    h = (length(v) - 1) * p + 1          # posição (pode ficar entre dois valores)
    i = floor(Int, h)
    i >= length(v) && return v[end]
    return v[i] + (h - i) * (v[i + 1] - v[i])
end

"""
    faixa_pam(curvas; inferior = 0.05, superior = 0.95)

Faixa da PAM a cada instante. `curvas` é uma lista de curvas de PAM (todas com
os mesmos instantes de tempo). Devolve `(inf, mediana, sup)`: em cada instante,
o percentil `inferior`, a mediana e o percentil `superior` entre as curvas.
Com o padrão (5% e 95%), a faixa contém 90% das curvas.
"""
function faixa_pam(curvas; inferior = 0.05, superior = 0.95)
    n = length(curvas[1])
    # valores de todas as curvas no instante k
    no_instante(k) = [c[k] for c in curvas]
    inf = [percentil(no_instante(k), inferior) for k in 1:n]
    med = [percentil(no_instante(k), 0.5) for k in 1:n]
    sup = [percentil(no_instante(k), superior) for k in 1:n]
    return (inf = inf, mediana = med, sup = sup)
end

"""
    periodos_de_risco(t, inf, sup; limiares = LIMIARES)

Divide o tempo em períodos. Em cada período, a borda de baixo da faixa (`inf`,
o pior caso) fica num mesmo nível de risco e a borda de cima (`sup`, o melhor
caso) fica num mesmo nível. Devolve uma lista de períodos com: `inicio` e `fim`
(s), `pior` e `melhor` (níveis) e `pam_min` e `pam_max` (mmHg) no período.
"""
function periodos_de_risco(t, inf, sup; limiares = LIMIARES)
    periodos = []
    for k in eachindex(t)
        pior, melhor = nivel(inf[k]; limiares), nivel(sup[k]; limiares)
        if !isempty(periodos) && periodos[end].pior == pior && periodos[end].melhor == melhor
            # mesmo par de níveis do instante anterior: estende o período atual
            p = periodos[end]
            periodos[end] = merge(p, (fim = t[k], pam_min = min(p.pam_min, inf[k]), pam_max = max(p.pam_max, sup[k])))
        else
            # mudou de nível: começa um período novo
            push!(periodos, (inicio = t[k], fim = t[k], pior = pior, melhor = melhor,
                             pam_min = inf[k], pam_max = sup[k]))
        end
    end
    return periodos
end

"""
    descrever(periodos)

Uma frase por período, por exemplo:
"de 45 a 600 s: PAM entre 58 e 80 mmHg — risco entre moderado e baixo".
"""
function descrever(periodos)
    map(periodos) do p
        risco = p.pior == p.melhor ? "risco $(p.pior)" : "risco entre $(p.pior) e $(p.melhor)"
        "de $(round(Int, p.inicio)) a $(round(Int, p.fim)) s: PAM entre " *
            "$(round(Int, p.pam_min)) e $(round(Int, p.pam_max)) mmHg — $risco"
    end
end
