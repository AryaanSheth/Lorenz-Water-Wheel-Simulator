main: main.cpp
	g++ -o main main.cpp -std=c++11 -ltrapfpe -lpgplot -lcpgplot -lX11 -lm
