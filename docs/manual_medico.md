# Manual da plataforma — previsão da PAM na indução com propofol e remifentanil

A plataforma monta um **Windkessel de 4 elementos** com os dados do paciente (a
árvore arterial dele, em analogia com um circuito elétrico) e simula, batimento a
batimento, a pressão arterial nos 20 minutos após o início do propofol. Mostra a
**faixa** em que a pressão arterial média (PAM) deve ficar, a PAS e a PAD,
**em que períodos** a PAM pode entrar em risco de hipotensão e quando o paciente
deve dormir (BIS abaixo de 60).

> Ferramenta de ensino e pesquisa. Não substitui a avaliação e a decisão
> clínica do anestesiologista.

## 1. Instalação (uma vez só)

1. **Instale o Julia.** Acesse [julialang.org/install](https://julialang.org/install/)
   e siga as instruções para Windows (o instalador pelo Microsoft Store é o mais
   simples).
2. **Baixe a plataforma.** Na página do projeto no GitHub, clique no botão verde
   **Code** e depois em **Download ZIP**.
3. **Descompacte** o ZIP numa pasta, por exemplo em Documentos.

## 2. Abrir a plataforma

1. Na pasta descompactada, dê **dois cliques** em **`abrir_plataforma.bat`**.
2. Abre uma janela preta (é o Julia trabalhando) e, em seguida, a plataforma no
   navegador. **Na primeira vez pode levar alguns minutos**; nas próximas, menos
   de um minuto.
3. **Não feche a janela preta** enquanto estiver usando. Para fechar a
   plataforma, feche a janela preta.

Se o navegador não abrir sozinho, digite no navegador: `http://127.0.0.1:8080`

## 3. Usar

1. **Paciente:** idade, sexo, peso e altura.
2. **Monitor antes da indução:** PAS, PAD, PAM e FC medidas antes de qualquer
   fármaco (antes do remifentanil). O **débito cardíaco** é opcional: se não for
   informado, é estimado pela superfície corporal (3,0 L/min/m²).
3. **Propofol:** escolha uma das duas formas:
   - **Dose:** dose em mg/kg, tempo de injeção e, se houver, a manutenção em
     mg/kg/h;
   - **Concentração-alvo:** a concentração plasmática desejada em µg/ml (como
     numa bomba de TCI).
4. **Remifentanil:** a taxa da infusão em µg/kg/min (0 = sem remifentanil) e
   quantos minutos antes do propofol ela começou. Os valores iniciais (0,12
   µg/kg/min, 9 min antes) são as medianas de 21 induções do VitalDB.
5. Clique em **Simular**. O cálculo leva alguns segundos; os gráficos de baixo
   aparecem depois.

Aceita vírgula ou ponto nos números (2,5 ou 2.5).

## 4. Como ler o resultado

- **Resumo (no alto):** o pior caso — até onde a PAM pode cair e em que período.
- **Gráfico:** o fundo mostra os níveis de risco:
  - verde: **baixo** (PAM ≥ 65 mmHg);
  - amarelo: **moderado** (55 a 65 mmHg);
  - vermelho: **alto** (PAM < 55 mmHg).

  A **faixa azul clara** é onde a PAM fica em 9 de cada 10 pacientes virtuais;
  a **escura**, em metade deles; a **linha azul** é a mediana. A **faixa roxa**
  no alto fica mais forte quanto mais pacientes estão dormindo (BIS < 60); as
  linhas roxas pontilhadas marcam quando, em média, o paciente dorme e acorda.
  O tempo está em segundos após o início do propofol.
- **PAS e PAD (linhas tracejadas):** pressões **centrais** (na aorta), que é o
  que o Windkessel calcula. A PAD e a PAM são praticamente as mesmas do braço; a
  PAS central é menor que a do manguito, porque a onda de pressão se amplifica
  até a periferia.
- **Windkessel do paciente:** os parâmetros da árvore arterial do paciente (Rp,
  C, Rc, L e a constante de tempo Rp·C) e a onda de pressão e de fluxo aórticos
  em três batimentos: o basal, o do início do propofol e o da PAM mínima.
- **Débito cardíaco** e **Como os fármacos mudam o Windkessel:** quanto a
  resistência (Rp), o volume sistólico e a FC mudam ao longo do tempo, em
  porcentagem do basal.
- **Pressão por dose de propofol:** a PAM mínima para cada dose de propofol,
  com o remifentanil fixo no valor escolhido.
- **Pressão por dose de remifentanil:** a PAM mínima para cada taxa de
  remifentanil, com o propofol fixo no valor escolhido.
- **Tabela de períodos:** para cada intervalo de tempo, entre quais valores a PAM
  fica e o nível de risco, do pior caso (borda de baixo da faixa) ao melhor
  (borda de cima). Exemplo: "257 a 1077 s — 49 a 68 mmHg — alto a baixo".

## 5. Limitações

- O quanto o propofol e o remifentanil mudam a resistência, o volume sistólico
  e a FC vem de um estudo em **voluntários saudáveis** (Su et al. 2023) e ainda
  **não foi validado em pacientes**.
- No Windkessel, Rc, L e a forma da onda de fluxo são os de Stergiopulos (1999)
  para todos os pacientes, e a complacência (C) fica constante durante a
  anestesia.
- A faixa usa a variabilidade publicada nos artigos; a variação real entre
  pacientes provavelmente é **maior**.
- O modelo considera **o propofol e o remifentanil**: não inclui outros opioides,
  outros hipnóticos, gases, vasopressores, volemia nem a laringoscopia.
- No modelo de BIS usado (Eleveld 2018), o remifentanil **não muda o BIS**
  diretamente.
- Todas as escolhas do modelo e as fontes de cada valor estão em
  [`docs/decisoes.md`](decisoes.md).

## 6. Problemas comuns

| O que aparece | O que fazer |
|---|---|
| "O Julia nao foi encontrado neste computador." | Instale o Julia (passo 1.1) e dê dois cliques de novo. |
| A página diz "Sem resposta do servidor". | A janela preta foi fechada. Dê dois cliques em `abrir_plataforma.bat` de novo. |
| "Preencha o campo: …" ou "Valor inválido em …" | Corrija o campo indicado e clique em Simular. |
| A janela preta mostra "address already in use" (endereço em uso). | A plataforma já está aberta em outra janela preta. Use aquela, ou feche-a e abra de novo. |
| A janela preta mostra um erro e não abre o navegador. | Tire um print da janela e envie para o grupo. |
