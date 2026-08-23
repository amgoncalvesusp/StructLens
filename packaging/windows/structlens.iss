#define MyAppName "StructLens"
#define MyAppPublisher "Adriano Marques Gonçalves (UNIARA)"
#ifndef MyAppVersion
  #define MyAppVersion "0.1.1"
#endif

[Setup]
AppId={{C7B67C28-3A83-4E45-8E67-6F4C3A6AA0D6}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppPublisherURL=https://github.com/amgoncalvesusp/StructLens
DefaultDirName={localappdata}\Programs\StructLens
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
OutputDir=..\..\release
OutputBaseFilename=StructLens-Setup
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
UninstallDisplayIcon={app}\StructLens.exe

[Files]
Source: "..\..\dist\StructLens.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\..\dist\structlens-*.whl"; DestDir: "{app}\packages"; Flags: ignoreversion
Source: "Install-StructLens.ps1"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\..\README.md"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\..\LICENSE"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{autoprograms}\StructLens"; Filename: "{app}\StructLens.exe"; Parameters: "--help"
Name: "{autodesktop}\StructLens"; Filename: "{app}\StructLens.exe"; Parameters: "--help"; Tasks: desktopicon

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; GroupDescription: "Additional shortcuts:"; Flags: unchecked

[Run]
Filename: "{app}\StructLens.exe"; Parameters: "--help"; Description: "Verify the StructLens installation"; Flags: postinstall nowait skipifsilent
