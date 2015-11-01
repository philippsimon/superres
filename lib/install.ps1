Function Download-And-Expand-ZIPFile($url, $path)
{
    $url -match '[^\/]+$' | Out-Null;
    $file = "$path\" + $matches[0]

    Write-Output "Create empty directory $path"
    if (Test-Path "$path")
    {
        Remove-Item -recurse -force -confirm:$false "$path" | Out-Null
    }
    New-Item -Path "$path" -Type directory | Out-Null

    Write-Output "Download $url"
    (New-Object System.Net.WebClient).DownloadFile($url, $file)

    Write-Output "Expand $file"
    Expand-ZIPFile "$file" "$path"

    Remove-Item -recurse -force -confirm:$false "$file" | Out-Null
}

$BITS = 64

$url = "http://heanet.dl.sourceforge.net/project/enblend/enblend-enfuse/enblend-enfuse-4.1/enblend-enfuse-4.1.4-win$BITS.zip"
Download-And-Expand-ZIPFile $url "$pwd\enblend-enfuse"

$url = "http://www.centrostudiprogressofotografico.it/download/dcraw-9.26-ms-$BITS-bit.zip"
$file = Download-And-Expand-ZIPFile $url "$pwd\dcraw"
Rename-Item "$pwd\dcraw\dcraw-9.26-ms-$BITS-bit.exe" "$pwd\dcraw\dcraw.exe"