RES=test_results_rdex
TEST=01.run_all_examples
FIX=fix_R

all: $(RES)

$(RES): $(TEST) $(FIX) $(wildcard rdex_*.R)
	$(TEST)

clean:
	rm -f *.Rout g.c.test .RData test_results_rdex
	rm -r plots
