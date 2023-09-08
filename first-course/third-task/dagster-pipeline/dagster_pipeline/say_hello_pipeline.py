from dagster import job

@job
def say_hello():
    print(f"Hello world from Dagster!")

if __name__ == '__main__':
    say_hello()