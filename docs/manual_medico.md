# Manual da plataforma — previsão da PAM na indução com propofol

A plataforma mostra a **faixa** em que a pressão arterial média (PAM) do paciente
deve ficar nos 20 minutos após o início do propofol, e **em que períodos** ela
pode entrar em risco de hipotensão.

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
2. **Monitor antes da indução:** PAM e FC medidas antes do propofol.
3. **Propofol:** escolha uma das duas formas:
   - **Dose:** dose em mg/kg, tempo de injeção e, se houver, a manutenção em
     mg/kg/h;
   - **Concentração-alvo:** a concentração plasmática desejada em µg/ml (como
     numa bomba de TCI).
4. Clique em **Prever**. O cálculo leva alguns segundos.

Aceita vírgula ou ponto nos números (2,5 ou 2.5).

## 4. Como ler o resultado

- **Resumo (no alto):** o pior caso — até onde a PAM pode cair e em que período.
- **Gráfico:** o fundo mostra os níveis de risco:
  - verde: **baixo** (PAM ≥ 65 mmHg);
  - amarelo: **moderado** (55 a 65 mmHg);
  - vermelho: **alto** (PAM < 55 mmHg).

  A **faixa azul** é onde a PAM fica em 90% das versões simuladas do paciente;
  a **linha azul** é a mediana. O tempo está em segundos após o início do
  propofol.
- **Tabela de períodos:** para cada intervalo de tempo, entre quais valores a PAM
  fica e o nível de risco, do pior caso (borda de baixo da faixa) ao melhor
  (borda de cima). Exemplo: "257 a 1077 s — 49 a 68 mmHg — alto a baixo".

## 5. Limitações

- O efeito do propofol na pressão vem de um estudo em **voluntários saudáveis**
  (Su et al. 2023) e ainda **não foi validado em pacientes**.
- A faixa usa a variabilidade publicada nos artigos; a variação real entre
  pacientes provavelmente é **maior**.
- O modelo considera **só o propofol**: não inclui opioides, outros hipnóticos,
  gases, vasopressores, volemia nem a laringoscopia.
- Todas as escolhas do modelo e as fontes de cada valor estão em
  [`docs/decisoes.md`](decisoes.md).

## 6. Problemas comuns

| O que aparece | O que fazer |
|---|---|
| "O Julia nao foi encontrado neste computador." | Instale o Julia (passo 1.1) e dê dois cliques de novo. |
| A página diz "Sem resposta do servidor". | A janela preta foi fechada. Dê dois cliques em `abrir_plataforma.bat` de novo. |
| "Preencha o campo: …" ou "Valor inválido em …" | Corrija o campo indicado e clique em Prever. |
| A janela preta mostra "address already in use" (endereço em uso). | A plataforma já está aberta em outra janela preta. Use aquela, ou feche-a e abra de novo. |
| A janela preta mostra um erro e não abre o navegador. | Tire um print da janela e envie para o grupo. |
