En aquesta carpeta es troben els arxius del segon apartat de la pràctica.
El arxiu de programa 'cami_sol.f90' al ser compilat calcula els angles que descriu la posició del Sol respecte una posició concreta a la Terra durant un any. La posició a la Terra queda determinada per l'angle de latitud 'alpha' que es defineix a l'inici del programa. El temps està discretitzat en dies (j) i minuts(i). j=0, i=0 es correspon al moment en què la Terra està al periheli (3 de Gener); j=0, i=1440 són les 24h d'aquest mateix dia.

Els resultats del càlcul queden escrits a l'arxiu 'angles.txt' i a l'arxiu 'dist_sol.txt'. El primer arxiu té una primera columna per l'angle horitzontal i una segona columna per l'angle vertical, i el segon arxiu conté la distància entre el Sol i la posició a la Terra per cada dia.

El programa també genera quatre arxius 'equinoci#.txt' amb les dades dels angles als equinocis i als solsticis. A partir d'aquests, al compilar el programa 'plot.gp' es genera un gràfic de la posició del Sol en l'arxiu 'equinocis.png'.
