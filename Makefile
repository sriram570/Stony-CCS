
all: poav2
	 
poav2:
	cd external/poaV2 && make poa 2>/dev/null

clean:
	cd external/poaV2 && make clean
