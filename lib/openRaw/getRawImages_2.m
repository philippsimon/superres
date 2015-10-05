function theRawImage = getRawImages_2( theFilename)

myWD = pwd();
[myPath, myName, myExt, myVer] = fileparts( theFilename);
myImageName = [ myName, myExt];
if strcmp( myPath, myWD)==0
    %fremdes Verzeichnis:
    copyfile( 'dcrawx86', myPath);
    fDel = 1;
else
    fDel = 0;
end

cd( myPath);

theRawImage = openRaw_2( myImageName);

if fDel
    delete( 'dcrawx86');
end

cd( myWD);
