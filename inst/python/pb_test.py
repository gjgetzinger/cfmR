def pb_test():
  import time 
  from tqdm import tqdm
  for i in tqdm(range(10)):
    time.sleep(3)
