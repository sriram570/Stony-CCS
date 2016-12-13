
all: poav2
	 
poav2:
	cd external/poav2 && make poa 2>/dev/null

clean:
	cd external/poav2 && make clean
