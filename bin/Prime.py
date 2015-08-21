from math import sqrt

class Prime:
    
    primes = []

    @staticmethod
    def init(maxPrime=100000):
        vals = [0]*maxPrime

        

        Prime.primes = []

        for i in range(2,maxPrime):

            maxDiv = int(sqrt(i))

            prime=True
            for j in range(2,maxDiv+1):
                if i%j==0:
                    prime=False
                    break

            if prime:
                Prime.primes.append(i)

    @staticmethod
    def getPrime(N):
        return Prime.primes[ N % len(Prime.primes) ]

if __name__ == '__main__':
    
    Prime.init(10)
#    print Prime.primes

    for i in range(20):
        print Prime.getPrime(i)
