#define MyAppName "StructLens"
#define MyAppPublisher "Adriano Marques Gonçalves (UNIARA)"
#ifndef MyAppVersion
#define MyAppVersion "0.3.0"
#endif

[Setup]
AppId={{C7B67C28-3A83-4E45-8E67-6F4C3A6AA0D6}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppPublisherURL=https://github.com/amgoncalvesusp/StructLens
SetupIconFile=structlens.ico
DefaultDirName={autopf}\StructLens
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
PrivilegesRequired=admin
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
OutputDir=..\..\release
OutputBaseFilename=StructLens-v{#MyAppVersion}-Windows-x86_64-Setup
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
UninstallDisplayIcon={app}\gui\StructLens.exe

[Files]
Source: "..\..\dist\StructLens\*"; DestDir: "{app}\gui"; Flags: recursesubdirs createallsubdirs ignoreversion
Source: "..\..\dist\structlens.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\..\dist\structlens-*.whl"; DestDir: "{app}\packages"; Flags: ignoreversion
Source: "Install-StructLens.ps1"; DestDir: "{app}"; Flags: ignoreversion
Source: "structlens.ico"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\..\README.md"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\..\LICENSE"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{autoprograms}\StructLens"; Filename: "{app}\gui\StructLens.exe"
Name: "{autoprograms}\StructLens CLI"; Filename: "{app}\structlens.exe"; Parameters: "--help"
Name: "{autodesktop}\StructLens"; Filename: "{app}\gui\StructLens.exe"; Tasks: desktopicon

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; GroupDescription: "Additional shortcuts:"; Flags: unchecked

[Run]
Filename: "{app}\gui\StructLens.exe"; Description: "Launch StructLens"; Flags: postinstall nowait skipifsilent
