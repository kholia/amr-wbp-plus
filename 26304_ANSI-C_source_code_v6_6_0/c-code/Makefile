dummy:
	@# i686-w64-mingw32-gcc -Wall -s -o decoder-new.exe common/*.c lib_amr/*.c decoder/*.c 3gplib/Release/er-libisomedia.dll -L3gplib/Release/ -l3gplib -lm
	@# i686-w64-mingw32-gcc -O2 -s -o decoder-new.exe common/*.c lib_amr/*.c decoder/*.c 3gplib/Release/er-libisomedia.dll -L3gplib/Release/ -l3gplib -lm
	gcc -Wno-incompatible-pointer-types -Wno-implicit-function-declaration -Wno-unused-result -Wno-multichar -O2 -s -o decoder-new common/*.c lib_amr/*.c decoder/*.c libbase/*.c -lm
	@# i686-w64-mingw32-gcc -s -O2 -o encoder-new.exe common/*.c lib_amr/*.c encoder/*.c 3gplib/Release/er-libisomedia.dll -L3gplib/Release/ -l3gplib -lm
