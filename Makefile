all:
	javac -cp .:jbotsim.jar *.java

%.dot: all
	java -cp .:jbotsim.jar Main $@

clean:
	rm *.class
