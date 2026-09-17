param(
	[ValidateSet("Release", "Debug")]
	[string]$Configuration = "Release",

	[switch]$ContinueOnFailure
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $PSScriptRoot
$exePath = Join-Path $repoRoot "build\stage\bin\$Configuration\udf_benchmark.exe"

$configs = @(
	"cgogn\geometry\apps\configs\udf_benchmark_batch_deepfashion_medium.json",
	"cgogn\geometry\apps\configs\udf_benchmark_batch_deepfashion_dense.json"
)

function Invoke-CheckedCommand {
	param(
		[Parameter(Mandatory = $true)]
		[string]$StepName,

		[Parameter(Mandatory = $true)]
		[scriptblock]$Command
	)

	Write-Host ""
	Write-Host "===== $StepName ====="
	& $Command
	if ($LASTEXITCODE -ne 0) {
		throw "$StepName failed with exit code $LASTEXITCODE"
	}
}

if (-not (Test-Path -LiteralPath $exePath)) {
	throw "udf_benchmark executable not found: $exePath"
}

$failedConfigs = New-Object System.Collections.Generic.List[string]

foreach ($relativeConfig in $configs) {
	$configPath = Join-Path $repoRoot $relativeConfig
	if (-not (Test-Path -LiteralPath $configPath)) {
		throw "Benchmark config not found: $configPath"
	}

	$stepName = "Run $relativeConfig"
	try {
		Invoke-CheckedCommand -StepName $stepName -Command {
			& $exePath $configPath
		}
	}
	catch {
		$failedConfigs.Add($relativeConfig)
		Write-Warning $_
		if (-not $ContinueOnFailure) {
			throw
		}
	}
}

if ($failedConfigs.Count -gt 0) {
	throw "Benchmark queue finished with failed configs: $($failedConfigs -join ', ')"
}

Write-Host ""
Write-Host "Benchmark queue finished successfully."
