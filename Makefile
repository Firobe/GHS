TESTS=tutte \
	  cycle\ 100 \
	  hypercube\ 5 \
	  planar\ 75\ 3\ 3 \

all:
	@javac -cp .:jbotsim.jar *.java

%.dot: all
	@java -cp .:jbotsim.jar Main $@

clean:
	rm *.class

gengraph:
	@gcc -O3 -o gengraph gengraph.c -lm

gentest: gengraph
	@./gengraph $(GRAPH) -format dot | dot -Tdot | sed 's/^\s//g' | sed 's/ //g' | sed 's/\s/ /g' > test.dot

test: gentest all test.dot

tests: $(TESTS)

$(TESTS):
	@./gengraph $@ -format dot | dot -Tdot | sed 's/^\s//g' | sed 's/ //g' | sed 's/\s/ /g' > test.dot && java -cp .:jbotsim.jar Main test.dot
