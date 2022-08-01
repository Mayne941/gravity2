import time 

def benchmark_start(name):
	print("*"*50)
	print(f"BENCHMARKING {name}")
	print("*"*50)
	return time.time()

def benchmark_end(name, start):
	print("*"*50)
	print(f"Benchmarking COMPLETE {name}")
	print(f"TIME TO COMPLETE: {time.time() - start}")
	print("*"*50)