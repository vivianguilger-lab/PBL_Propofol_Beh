# PropofolWK

Simulador da **pressão arterial** e da **profundidade anestésica (BIS)** na indução da anestesia com
**propofol** e, se o médico escolher, **remifentanil**. A previsão é mostrada como **faixas**, que
representam a variabilidade entre pacientes, com os **períodos de risco de hipotensão** e gráficos de
**pressão por dose**, em uma plataforma local que o médico usa sem programar.

![Plataforma: ficha do paciente, nível de risco e indicadores](docs/imagens/plataforma.png)

> Ferramenta de ensino e pesquisa (PBL, 6º semestre de Engenharia Biomédica). Os efeitos na
> circulação vêm de um estudo em voluntários saudáveis e o modelo ainda não foi comparado à pressão
> medida em pacientes: **não use para decisões clínicas.**

## Sumário

- [O que o simulador faz](#o-que-o-simulador-faz)
- [Como instalar e abrir](#como-instalar-e-abrir)
- [O que a plataforma mostra](#o-que-a-plataforma-mostra)
- [Como o modelo funciona](#como-o-modelo-funciona)
- [Dados do VitalDB](#dados-do-vitaldb)
- [Verificação](#verificação)
- [Exemplo: pacientes A, B e C](#exemplo-pacientes-a-b-e-c)
- [Uso no Julia e estrutura do código](#uso-no-julia-e-estrutura-do-código)
- [Limitações e próximos passos](#limitações-e-próximos-passos)
- [Documentação](#documentação)
- [Referências](#referências)
- [Equipe](#equipe)

## O que o simulador faz

A partir dos dados do paciente e do protocolo escolhido, calcula segundo a segundo, nos 20 minutos
após o início do propofol:

- a **pressão arterial média (PAM)** e a **frequência cardíaca (FC)**;
- o **BIS** e quando o paciente dorme (BIS < 60) e acorda;
- as **concentrações** de propofol e de remifentanil no sangue;
- os **períodos de risco**: baixo (PAM ≥ 65 mmHg), moderado (55 a 64) e alto (< 55);
- a **PAM mínima para cada dose** de propofol e de remifentanil.

Tudo é calculado para **500 pacientes virtuais**: versões do mesmo paciente com os parâmetros sorteados
dentro da variabilidade publicada nos artigos. Só entram modelos e valores publicados e verificados; cada
escolha está registrada, com a fonte, em [`docs/decisoes.md`](docs/decisoes.md).

## Como instalar e abrir

1. Instale o **Julia** (1.10 ou mais novo): [julialang.org/install](https://julialang.org/install/).
2. Baixe o projeto: botão verde **Code → Download ZIP** (ou `git clone`) e descompacte.
3. Dê **dois cliques** em **`abrir_plataforma.bat`** (Windows). Abre uma janela preta (o Julia) e, em
   seguida, a plataforma no navegador. Na primeira vez pode levar alguns minutos.
4. Não feche a janela preta enquanto usa a plataforma. Se o navegador não abrir, acesse
   `http://127.0.0.1:8080`.

Tudo roda no próprio computador: nenhum dado sai da máquina. O pacote usa só bibliotecas que já vêm
com o Julia (LinearAlgebra, Random e Sockets), sem nada para baixar da internet.

## O que a plataforma mostra

### Ficha do paciente

| Seção | Campos |
|---|---|
| 1. Paciente | Idade, sexo, peso, altura |
| 2. Monitor | PAM e FC **antes de qualquer fármaco** (antes do remifentanil) |
| 3. Protocolo de indução | **Dose** (mg/kg, tempo de injeção e manutenção em µg/kg/min) ou **concentração-alvo** no plasma (µg/ml, como numa bomba de TCI) |
| 4. Remifentanil | Taxa (µg/kg/min; 0 = sem remifentanil) e quantos minutos antes do propofol começou. Valores iniciais: 0,12 µg/kg/min e 9 min, as medianas de 21 induções do VitalDB |

<img src="docs/imagens/ficha_medica.png" alt="Ficha médica preenchida" width="320">

### Resultados

No alto: o **nível de risco do pior caso** (5% dos pacientes virtuais ficam abaixo), uma régua com a PAM
basal, a mínima da mediana e a do pior caso, e seis indicadores (PAM basal → mínima, queda, sono, tempo com
PAM < 65 e < 55, concentração máxima). Abaixo, oito cartões:

**01 · Pressão arterial e sono (gráfico principal).** Mediana da PAM, faixa de metade dos pacientes,
faixa de 9 em cada 10, zonas de risco, a faixa roxa no alto (quanto mais forte, mais pacientes dormindo)
e quando, em média, o paciente dorme e acorda. Embaixo, a tabela de períodos de risco.

![Pressão arterial e sono](docs/imagens/pressao_e_sono.png)

**02 · Pressão por dose de propofol** (remifentanil fixo) e **03 · Pressão por dose de remifentanil**
(propofol fixo). PAM mínima do paciente típico para cada dose, sobre as zonas de risco, com as doses em que
o risco muda e a dose a partir da qual o paciente, em média, dorme. Clicar no gráfico simula aquela dose.

![Pressão por dose de propofol](docs/imagens/pressao_por_dose_propofol.png)

![Pressão por dose de remifentanil](docs/imagens/pressao_por_dose_remifentanil.png)

**04 · Profundidade anestésica (BIS).** Mediana, faixas, faixa-alvo de 40 a 60 e quando, em média, o
paciente dorme e acorda. O BIS aparece como o monitor mostra, com o atraso dele.

![BIS](docs/imagens/bis.png)

**05 · Frequência cardíaca**, **06 · Concentrações no plasma** (propofol e remifentanil),
**07 · Efeitos dos fármacos** (resistência periférica, volume sistólico e FC em % do basal) e
**08 · Paciente virtual** (todos os parâmetros usados, com a fonte):

![Paciente virtual](docs/imagens/paciente_virtual.png)

No rodapé: imprimir/salvar PDF e baixar as curvas em CSV. Há tema claro e escuro, e a página funciona
no celular.

<img src="docs/imagens/celular.png" alt="Plataforma no celular" width="300">

## Como o modelo funciona

```
propofol (mg/s)     → Eleveld 2018 → Cp (µg/ml) ─┬→ Su 2023 → PAM, FC, RPT, VS
                                                 └→ sítio de efeito → Ce → BIS (Eleveld 2018)
remifentanil (µg/s) → Minto 1997  → Cr (ng/ml) ──→ Su 2023
```

Passo de 1 s; t = 0 é o início do propofol, e o remifentanil começa antes (a conta começa em t < 0).

| Parte | Modelo | Fonte |
|---|---|---|
| Farmacocinética do propofol | 3 compartimentos; termos de opioide quando há remifentanil | Eleveld et al., 2018 |
| Farmacocinética do remifentanil | 3 compartimentos, massa magra de James | Minto et al., 1997 |
| Circulação | Renovação de RPT, VS e FC com feedback da PAM (barorreflexo) e interação propofol–remifentanil | Su et al., 2023 |
| Sono | Sítio de efeito (ke0), curva de Hill e atraso do monitor | Eleveld et al., 2018 |
| Variabilidade | 500 pacientes virtuais com as variâncias publicadas | Mesmos artigos |
| Risco | 65 e 55 mmHg | Salmasi et al., 2017; Walsh et al., 2013 |

**Farmacocinética.** Três compartimentos (sangue, rápido e lento), resolvidos de forma exata a cada segundo:

```math
\frac{dC_1}{dt} = -(k_{10}+k_{12}+k_{13})\,C_1 + k_{12}C_2 + k_{13}C_3 + \frac{u(t)}{V_1}, \qquad \frac{dC_i}{dt} = k_{i1}\,(C_1 - C_i)
```

**BIS.** O cérebro segue o sangue com atraso, e o BIS cai com a concentração no cérebro:

```math
\frac{dC_e}{dt} = k_{e0}\,(C_p - C_e), \qquad BIS = 93\left(1 - \frac{C_e^{\gamma}}{C_{e50}^{\gamma} + C_e^{\gamma}}\right)
```

com Ce50 = 3,08·e^(−0,00635·(idade − 35)) µg/ml, γ = 1,89 abaixo do Ce50 e 1,47 acima, e o atraso do
monitor de 15 + e^(0,0517·idade) s.

**Circulação.** PAM = RPT × VS × FC, cada um como fração do basal:

```math
\frac{dX}{dt} = k_{out}\left[\left(\frac{PAM}{PAM_0}\right)^{fb}\left(1 + E^{P}_{X}\right) - X\left(1 - E^{R}_{X}\right)\right], \qquad X \in \{RPT,\ VS,\ FC\}
```

k_out = 0,072/min e fb = −0,661. E^P é o efeito do propofol (vasodilatação e queda do volume sistólico,
maior nos idosos) e E^R, o do remifentanil, com os termos de interação do artigo. Integração por
Runge-Kutta de 4ª ordem.

**Pacientes virtuais e faixas.** Cada parâmetro com variabilidade publicada é sorteado (log-normal:
valor típico × e^η, η ~ Normal(0, ω²); três do remifentanil são aditivos, como no PAS). A cada segundo,
os percentis 5–95% e 25–75% das 500 curvas formam as faixas.

**Risco.** O pior caso (percentil 5%), com a PAM arredondada para mmHg inteiros, como no monitor, é
classificado segundo a segundo; o nível do paciente é o do pior segundo.

Explicação completa, com todas as equações, valores e fontes:
[`docs/como_o_modelo_funciona.html`](docs/como_o_modelo_funciona.html) (baixe e abra no navegador).

## Dados do VitalDB

O [VitalDB](https://vitaldb.net) é um banco aberto com 6.388 cirurgias. A seleção de casos (21 induções
com pressão invasiva e bomba de propofol) foi feita pelo grupo no repositório
[PBL_6](https://github.com/vivianguilger-lab/PBL_6). Aqui eles foram usados para:

- **descrever o remifentanil na prática:** em todas as 21 induções ele começou antes do propofol (mediana de
  9 min antes), com taxa média de 0,12 µg/kg/min. São os valores iniciais da plataforma;
- **conferir a farmacocinética de Minto:** com a taxa registrada pela bomba (`Orchestra/RFTN20_RATE`), a
  concentração de Minto diferiu da que a própria bomba mostrou (`Orchestra/RFTN20_CP`) em **0,97%**
  (mediana).

Tabela por caso: [`dados/remifentanil_vitaldb.csv`](dados/remifentanil_vitaldb.csv); script:
[`scripts/remifentanil_vitaldb.py`](scripts/remifentanil_vitaldb.py).

## Verificação

As equações foram conferidas contra o [Python Anesthesia Simulator](https://github.com/AnesthesiaSimulation/Python_Anesthesia_Simulator)
(PAS, commit 4afc03a), em 109 testes automáticos:

| O que | Diferença para o PAS |
|---|---|
| Concentração de propofol (com e sem opioide) e de remifentanil | < 10⁻¹³ |
| PAM e FC, só propofol e com os dois fármacos | < 10⁻¹² mmHg e bpm |
| Concentração no cérebro e equação do BIS | idênticas (< 1 ponto de BIS com o atraso exato do artigo) |

Os testes também conferem a concentração-alvo, o sorteio dos pacientes virtuais (média e variância), as
faixas e os períodos de risco e a plataforma. As tabelas de referência são geradas pelos scripts
`test/gerar_referencia_*.py`.

## Exemplo: pacientes A, B e C

Propofol em bolus de 30 s seguido de manutenção de 100 µg/kg/min, com e sem remifentanil
(0,12 µg/kg/min desde 9 min antes):

| Paciente | Remifentanil | PAM mínima, mediana (mmHg) | Pior caso (mmHg) | Risco | Dorme aos |
|---|---|---|---|---|---|
| A: mulher, 62 anos, 72 kg, PAM 104, 1,5 mg/kg | Sem | 77 (−26%) | 60 | Moderado | 149 s |
| | Com | 58 (−44%) | 47 | Alto | 131 s |
| B: homem, 28 anos, 80 kg, PAM 89, 2,0 mg/kg | Sem | 72 (−19%) | 57 | Moderado | 120 s |
| | Com | 58 (−35%) | 47 | Alto | 112 s |
| C: homem, 81 anos, 64 kg, PAM 107, 1,5 mg/kg | Sem | 66 (−38%) | 51 | Alto | 167 s |
| | Com | 48 (−55%) | 39 | Alto | 155 s |

Os três estão nos botões de exemplo da plataforma.

## Uso no Julia e estrutura do código

```julia
] activate .
] test                      # roda os testes

using PropofolWK
r = prever(idade = 62, peso = 72, altura = 160, sexo = "F", pam = 104, fc = 78,
           dose_mgkg = 1.5, duracao_bolus_s = 30, manutencao_mgkgh = 6,
           remi_ugkgmin = 0.12, remi_antes_s = 540)
r.metricas                  # PAM mínima, pior caso, risco, tempo abaixo de 65 e 55...
r.texto                     # períodos de risco em frases
plataforma()                # abre a plataforma
```

| Arquivo | O que faz |
|---|---|
| `src/farmacocinetica.jl` | Eleveld 2018: parâmetros do paciente e solução exata dos 3 compartimentos |
| `src/alvo.jl` | Concentração-alvo (TCI no plasma) |
| `src/remifentanil.jl` | Minto 1997 |
| `src/bis.jl` | Sítio de efeito e BIS de Eleveld |
| `src/hemodinamica.jl` | Modelo de Su (propofol, remifentanil e interações) |
| `src/populacao.jl` | Pacientes virtuais |
| `src/risco.jl` | Faixas, níveis e períodos de risco |
| `src/previsao.jl` | `prever`, pressão por dose (`risco_por_dose`, `risco_por_remi`) |
| `src/servidor.jl` | Servidor local da plataforma |
| `src/grafico.jl` | Gráfico da faixa em SVG |
| `app/index.html` | A página da plataforma |
| `test/` | Testes e scripts das tabelas de referência |
| `exemplos/prever_paciente.jl` | Exemplo de uso |

## Limitações e próximos passos

- Os efeitos na circulação foram estimados em voluntários saudáveis; falta comparar com a PAM medida nos
  casos do VitalDB e calibrar os parâmetros mais influentes.
- Não entram laringoscopia, estímulo cirúrgico, fluidos, vasopressores, outros opioides e hipnóticos.
- Só parte dos parâmetros tem variabilidade publicada: as faixas provavelmente são mais estreitas que a realidade.
- No modelo de BIS de Eleveld, o remifentanil não entra diretamente no BIS.
- O risco é definido pelo pior segundo, sem tempo mínimo abaixo do limiar.

Pendências detalhadas em [`docs/decisoes.md`](docs/decisoes.md).

## Documentação

- [`docs/como_o_modelo_funciona.html`](docs/como_o_modelo_funciona.html): tudo sobre o modelo: matemática,
  valores, fontes, saídas, cada parte da plataforma e o código.
- [`docs/manual_medico.md`](docs/manual_medico.md): manual da plataforma para o médico.
- [`docs/decisoes.md`](docs/decisoes.md): cada decisão de método, com a fonte e a situação.

## Referências

- Eleveld DJ et al. Pharmacokinetic–pharmacodynamic model for propofol for broad application in anaesthesia and sedation. *Br J Anaesth* 2018;120(5):942–59. doi:10.1016/j.bja.2018.01.018
- Su H et al. Pharmacodynamic mechanism-based interaction model for the haemodynamic effects of remifentanil and propofol in healthy volunteers. *Br J Anaesth* 2023;131(2):222–33. doi:10.1016/j.bja.2023.04.043
- Minto CF et al. Influence of age and gender on the pharmacokinetics and pharmacodynamics of remifentanil: I. Model development. *Anesthesiology* 1997;86(1):10–23.
- Salmasi V et al. Relationship between intraoperative hypotension, defined by either reduction from baseline or absolute thresholds, and acute kidney and myocardial injury after noncardiac surgery. *Anesthesiology* 2017;126(1):47–65.
- Walsh M et al. Relationship between intraoperative mean arterial pressure and clinical outcomes after noncardiac surgery. *Anesthesiology* 2013;119(3):507–15.
- NICE. Depth of anaesthesia monitors – Bispectral Index (BIS), E-Entropy and Narcotrend-Compact M. Diagnostics guidance DG6, 2012.
- Lee HC et al. VitalDB, a high-fidelity multi-parameter vital signs database in surgical patients. *Sci Data* 2022;9:279.
- Aubouin-Pairault B et al. PAS: a Python Anesthesia Simulator for drug control. *J Open Source Softw* 2023;8(88):5480.

## Equipe

Grupo 5 · Engenharia Biomédica · Faculdade Israelita de Ciências da Saúde Albert Einstein · Project Based
Learning, 6º semestre (2026). Orientador: Federico Aletti.

Ana Paula Villela Segura, Beatriz Neves de Oliveira, Gustavo Frug Mauro, Leticia Kairuz, Maria Gabriela
Medeiros von Adamek e Vivian Polakiewicz Guilger.

Licença: GPL-3.0 (ver [`LICENSE`](LICENSE)).
