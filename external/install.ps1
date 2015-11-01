function Expand-ZIPFile($file, $destination)
{
    $shell = new-object -com shell.application
    $zip = $shell.NameSpace($file)
    foreach($item in $zip.items())
    {
        $shell.Namespace($destination).copyhere($item)
    }
}

Function Download-And-Expand-ZIPFile($url, $path)
{
    $url -match '[^\/]+$' | Out-Null;
    $file = "$path\" + $matches[0]

    Write-Debug "Create empty directory $path"
    if (Test-Path "$path")
    {
        Remove-Item -recurse -force -confirm:$false "$path" | Out-Null
    }
    New-Item -Path "$path" -Type directory | Out-Null

    Write-Debug "Download $url"
    (New-Object System.Net.WebClient).DownloadFile($url, $file)

    Write-Debug "Expand $file"
    Expand-ZIPFile "$file" "$path"

    Remove-Item -recurse -force -confirm:$false "$file" | Out-Null
}

$ENV_BITS = [IntPtr]::Size * 8

Write-Output "# Install dcraw"
$url = "http://www.centrostudiprogressofotografico.it/download/dcraw-9.26-ms-$ENV_BITS-bit.zip"
$file = Download-And-Expand-ZIPFile $url "$pwd\dcraw"
Rename-Item "$pwd\dcraw\dcraw-9.26-ms-$ENV_BITS-bit.exe" "$pwd\dcraw\dcraw.exe"

Write-Output "# Install enblend-enfuse"
$url = "http://heanet.dl.sourceforge.net/project/enblend/enblend-enfuse/enblend-enfuse-4.1/enblend-enfuse-4.1.4-win$ENV_BITS.zip"
Download-And-Expand-ZIPFile $url "$pwd\enblend-enfuse"

Write-Output "# Install tufuse"
if  ($ENV_BITS -eq 64)
{
    $url = "http://www.tawbaware.com/tufusex64.zip"
}
else
{
    $url = "http://www.tawbaware.com/tufuse.zip"
}
$file = Download-And-Expand-ZIPFile $url "$pwd\tufuse"