### How to build（Windows）

cd imagej-plugins\src
build.bat C:\path\to\ij.jar

Then `SP8CT_Plugins.jar` is generated.


### Imstall

1. Put `SP8CT_Plugins.jar` in `plugins/` folder.

### How to open with D&D

** For ImageJ. **  
Extract `HandleExtraFileTypes.class` from the `build/` directory and place it
directly in the root of the `plugins/` folder (as a standalone `.class` file,
not inside a JAR). ImageJ will automatically call this class for unknown file extensions.

** For Fiji. **  
Since Fiji already has its own `HandleExtraFileTypes`, the method described above cannot be used.
Instead, please select a file from `Plugins > SP8CT > IO > Open HIS/IMG...`.

