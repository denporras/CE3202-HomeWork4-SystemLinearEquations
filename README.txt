Como Framework de pruebas unitarias se empleo GTest de Google.

Instrucciones:

Instalar GTest: Seguir los siguientes comandos

	sudo apt-get install libgtest-dev
	sudo apt-get install cmake # install cmake
	cd /usr/src/gtest
	sudo cmake CMakeLists.txt
	sudo make
	sudo cp *.a /usr/lib

Para compilar el programa con cmake, correr en terminal:
	cmake CMakeLists.txt
	make

Para correr las pruebas unitarias ejecutar en terminal:
	./runTests
Para correr el programa y probar las funciones ejecutar en terminal: 
	./tarea
	
Seguir las instrucciones que indica el programa, las matrices correspondiente a cada numero estan acontinuacion.

*******************IMPORTANTE, MATRICES CON SU NUMERO Y DESCRIPCION******************
Matrices:

1. 	|4,5,7|     	MATRIZ NORMAL
  	|1,2,4|
   	|1,1,3|

2.	|6, 1,9|    	MATRIZ NORMAL
	|8,14,3|
	|3, 2,7|

3. 	|1   ,0.01|   	MATRIZ MAL CONDICIONADA
	|0.99,1   |

4. 	|4  ,5|  	MATRIZ MAL CONDICIONADA
	|4.1,5|

5. 	|1,1,1|  	MATRIZ NO INVERTIBLE  
	|1,1,1|
	|1,1,1|

Sistemas de ecuaciones:

1.	  |1,1,3|
	A=|1,2,4|
	  |4,5,7|

	b=|7,11,20|

2.
	  |1,0,5,7,0,7 | 
	  |3,3,44,6,8,8| 
	A=|6,4,3,3,4,6 | 
	  |89,8,7,6,5,4| 
	  |5,6,8,8,9,8 | 
	  |65,4,3,3,4,5| 

	b=|4,7,7,4,7,7|

3.
	  |2,1,3,5  | DIMENSIONES INCORRECTAS
	  |-1,0,7,1 |
	A=|0,-1,-1,3|
     	  |-3,7,4,3 |
	  |1,6,4,-3 |

	b=|7,11,20|

