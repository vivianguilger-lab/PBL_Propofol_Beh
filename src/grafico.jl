# =============================================================================
# grafico.jl — gráfico da faixa da PAM sobre as zonas de risco
#
# Gera um arquivo SVG (imagem que abre em qualquer navegador) sem pacote de
# gráficos: o desenho é escrito direto, elemento por elemento.
#   - fundo dividido nas zonas de risco (baixo, moderado, alto);
#   - faixa da PAM (área) e mediana (linha);
#   - eixos com "grandeza (unidade)", tempo em segundos, sem quadriculado.
# =============================================================================

# Cores das zonas de risco, da faixa e do texto.
const CORES = (baixo = "#e6f4ea", moderado = "#fff3dc", alto = "#fde7e7",
               faixa = "#2a78d6", texto = "#333333", eixo = "#888888")

# Leva um valor do intervalo [a, b] para o intervalo [p0, p1] (posição em pixels).
escala(v, a, b, p0, p1) = p0 + (v - a) / (b - a) * (p1 - p0)

# Escapa os caracteres especiais do SVG (&, < e >) num texto.
escapar(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")

"""
    svg_faixa(t, faixa; titulo = "", limiares = LIMIARES)

Gráfico da faixa da PAM ao longo do tempo, como texto SVG.
`t` é o tempo (s) e `faixa` vem de `faixa_pam` (`inf`, `mediana`, `sup`).
"""
function svg_faixa(t, faixa; titulo = "", limiares = LIMIARES)
    # tamanho da imagem e margens (pixels)
    L, A = 900, 500
    esq, dir, topo, base = 70, 110, 50, 60
    # limites dos eixos: tempo de 0 ao fim; PAM com folga, em múltiplos de 10
    t0, t1 = first(t), last(t)
    y0 = 10 * floor(min(minimum(faixa.inf), limiares.alto) / 10 - 0.5)
    y1 = 10 * ceil(maximum(faixa.sup) / 10 + 0.5)
    X(v) = escala(v, t0, t1, esq, L - dir)
    Y(v) = escala(v, y0, y1, A - base, topo)        # no SVG o y cresce para baixo

    s = IOBuffer()
    println(s, """<svg xmlns="http://www.w3.org/2000/svg" width="$L" height="$A" viewBox="0 0 $L $A" font-family"sans-serif" font-size="13">""")
    println(s, """<rect width="$L" height="$A" fill="white"/>""")

    # zonas de risco: retângulos de fundo, com o nome à direita
    zonas = [(:baixo, limiares.moderado, y1), (:moderado, limiares.alto, limiares.moderado), (:alto, y0, limiares.alto)]
    for (nome, ya, yb) in zonas
        println(s, """<rect x="$(X(t0))" y="$(Y(yb))" width="$(X(t1) - X(t0))" height="$(Y(ya) - Y(yb))" fill="$(getfield(CORES, nome))"/>""")
        println(s, """<text x="$(X(t1) + 8)" y="$((Y(ya) + Y(yb)) / 2 + 4)" fill="$(CORES.texto)">risco $nome</text>""")
    end

    # faixa: contorno pela borda de cima (ida) e pela de baixo (volta)
    ida = join(["$(X(t[k])),$(Y(faixa.sup[k]))" for k in eachindex(t)], " ")
    volta = join(["$(X(t[k])),$(Y(faixa.inf[k]))" for k in reverse(eachindex(t))], " ")
    println(s, """<polygon points="$ida $volta" fill="$(CORES.faixa)" fill-opacity="0.3" stroke="none"/>""")
    # mediana
    med = join(["$(X(t[k])),$(Y(faixa.mediana[k]))" for k in eachindex(t)], " ")
    println(s, """<polyline points="$med" fill="none" stroke="$(CORES.faixa)" stroke-width="2"/>""")

    # eixos: linhas, marcas e números (tempo a cada 120 s; PAM a cada 10 mmHg)
    println(s, """<line x1="$esq" y1="$(A - base)" x2="$(L - dir)" y2="$(A - base)" stroke="$(CORES.eixo)"/>""")
    println(s, """<line x1="$esq" y1="$topo" x2="$esq" y2="$(A - base)" stroke="$(CORES.eixo)"/>""")
    for v in 0:120:t1
        println(s, """<line x1="$(X(v))" y1="$(A - base)" x2="$(X(v))" y2="$(A - base + 5)" stroke="$(CORES.eixo)"/>""")
        println(s, """<text x="$(X(v))" y="$(A - base + 20)" text-anchor="middle" fill="$(CORES.texto)">$(round(Int, v))</text>""")
    end
    for v in y0:10:y1
        println(s, """<line x1="$(esq - 5)" y1="$(Y(v))" x2="$esq" y2="$(Y(v))" stroke="$(CORES.eixo)"/>""")
        println(s, """<text x="$(esq - 9)" y="$(Y(v) + 4)" text-anchor="end" fill="$(CORES.texto)">$(round(Int, v))</text>""")
    end
    # nomes dos eixos, no formato "grandeza (unidade)"
    println(s, """<text x="$((esq + L - dir) / 2)" y="$(A - 15)" text-anchor="middle" fill="$(CORES.texto)">Tempo após o início do bolus (s)</text>""")
    println(s, """<text x="20" y="$((topo + A - base) / 2)" text-anchor="middle" fill="$(CORES.texto)" transform="rotate(-90 20 $((topo + A - base) / 2))">PAM (mmHg)</text>""")

    # título e legenda
    println(s, """<text x="$esq" y="25" font-size="15" fill="$(CORES.texto)">$(escapar(titulo))</text>""")
    println(s, """<rect x="$(L - dir - 300)" y="$(topo + 10)" width="16" height="10" fill="$(CORES.faixa)" fill-opacity="0.3"/>""")
    println(s, """<text x="$(L - dir - 278)" y="$(topo + 19)" fill="$(CORES.texto)">faixa (90% das versões)</text>""")
    println(s, """<line x1="$(L - dir - 300)" y1="$(topo + 35)" x2="$(L - dir - 284)" y2="$(topo + 35)" stroke="$(CORES.faixa)" stroke-width="2"/>""")
    println(s, """<text x="$(L - dir - 278)" y="$(topo + 39)" fill="$(CORES.texto)">mediana</text>""")

    println(s, "</svg>")
    return String(take!(s))
end

# grafico_faixa: grava o gráfico em `arquivo` (.svg) e devolve o caminho do arquivo.
grafico_faixa(arquivo, t, faixa; kw...) = (write(arquivo, svg_faixa(t, faixa; kw...)); arquivo)
