all:
	gcc -o seedot seedot.c lib/libsike.a

clean:
	rm seedot

siot: 
	gcc -o siot postQuantumOT.c lib/libsike.a

slean: 
	rm siot

try:
	rm siot
	gcc -o siot postQuantumOT.c lib/libsike.a
	./siot
