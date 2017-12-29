all:
	javac -cp .:jbotsim.jar *.java

run: all
	java -cp .:jbotsim.jar HelloWorld

clean:
	rm *.class
