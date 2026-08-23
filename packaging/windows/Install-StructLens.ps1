[CmdletBinding()]
param(
    [string]$Python = "python"
)

$ErrorActionPreference = "Stop"
$bundleRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$wheel = Get-ChildItem -LiteralPath $bundleRoot -Filter "structlens-*.whl" | Select-Object -First 1

if ($null -eq $wheel) {
    throw "The StructLens wheel was not found beside this installer."
}

Write-Host "Installing StructLens from $($wheel.Name)..."
& $Python -m pip install --upgrade $wheel.FullName
if ($LASTEXITCODE -ne 0) {
    throw "pip could not install StructLens. Check the Python executable passed with -Python."
}

$pluginPath = & $Python -c "import pathlib, structlens; print(pathlib.Path(structlens.__file__).parent / 'plugin' / 'entrypoint.py')"
Write-Host "Installation complete."
Write-Host "PyMOL plugin entry point: $pluginPath"
Write-Host "In PyMOL, load that entry point through Plugin > Install Plugin."
