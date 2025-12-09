import subprocess

haveMPL=False   # do we have matplotlib?

try:
    #p = subprocess.call(['python3', '-c', 'import matplotlib;print(hasattr("matplotlib", "__version__"))'] )
    p = subprocess.call(['python3', '-c', 'import matplotlib;matplotlib.__version__'] )
    print(p)
    if p < 1:
       haveMPL=True
    #mplversion=p.stdout.readline().strip()
    #print('matplot version f{mplversion} found.')
    #haveMPL=True
except Exception as e:
    print(e)
    haveMPL=False

if not haveMPL:
    print("matplotlib not found. Testing of some examples will be skipped.")
