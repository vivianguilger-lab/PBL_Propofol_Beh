@echo off
rem ==========================================================================
rem abrir_plataforma.bat - abre a plataforma de previsao da PAM (PropofolWK)
rem De dois cliques neste arquivo. O navegador abre sozinho com a plataforma.
rem Para fechar a plataforma, feche esta janela preta.
rem (Sem acentos nas mensagens: o prompt do Windows pode nao mostra-los.)
rem ==========================================================================
title Plataforma - previsao da PAM com propofol

rem vai para a pasta onde este arquivo esta (a pasta do projeto)
cd /d "%~dp0"

rem confere se o Julia esta instalado
where julia >nul 2>nul
if errorlevel 1 (
    echo.
    echo O Julia nao foi encontrado neste computador.
    echo Instale pelo site https://julialang.org/install/ e de dois cliques aqui de novo.
    echo.
    pause
    exit /b 1
)

echo.
echo Abrindo a plataforma... Na primeira vez pode levar alguns minutos.
echo NAO feche esta janela enquanto estiver usando a plataforma.
echo.

rem prepara o projeto (resolve + instantiate) e abre a plataforma
julia --project=. -e "using Pkg; Pkg.resolve(); Pkg.instantiate(); using PropofolWK; plataforma()"

rem se o Julia parar com erro, a janela fica aberta para mostrar a mensagem
pause
