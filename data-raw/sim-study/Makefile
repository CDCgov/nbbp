rvec := 0.1 0.25 0.5 0.75 0.95 1.0 1.05 1.25 2.0
kvec := 0.1111111 0.3333333 1.0000000 3.0000000
nvec := 10 20 40
ivec := 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99

.SILENT: k.txt r.txt n.txt i.txt clean
.PHONY: simulate bayes maxlik check_timeout

all: bayes maxlik

maxlik: i.txt n.txt r.txt k.txt simulate
	src/maxlik.sh

bayes: i.txt n.txt r.txt k.txt simulate
	src/bayes.sh

# timeout is a downstream dependency but there's no point simulating if we can't analyze
simulate: i.txt n.txt r.txt k.txt check_timeout
	src/simulate.sh

i.txt:
	printf "%s\n" $(ivec) > par/i.txt

n.txt:
	printf "%s\n" $(nvec) > par/n.txt

r.txt:
	printf "%s\n" $(rvec) > par/r.txt

k.txt:
	printf "%s\n" $(kvec) > par/k.txt

# Makes sure timeout command exists, otherwise we can't use bayes.sh or maxlik.sh
check_timeout:
	@which timeout > which_timeout.txt

clean:
	rm which_timeout.txt
	rm -f i.txt
	rm -f n.txt
	rm -f r.txt
	rm -f k.txt
	rm -rf par
	rm -rf data
	rm -rf bayes
	rm -rf maxlik
	mkdir par
	mkdir data
	mkdir bayes
	mkdir maxlik
	cat <<EOF > par/.placeholder
	cat <<EOF > data/.placeholder
	cat <<EOF > bayes/.placeholder
	cat <<EOF > maxlik/.placeholder
