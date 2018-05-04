# -*- mode: python -*-

block_cipher = None


a = Analysis(['ScaffoldTree.py', 'ScaffoldTree.spec'],
             pathex=['/home/npapado/Documents/repos/merlot/inst/python'],
             binaries=[],
             datas=[],
             hiddenimports=['scipy._lib.messagestream',
                            'pandas._libs.tslibs.timedeltas'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='ScaffoldTree',
          debug=False,
          strip=False,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='ScaffoldTree')
