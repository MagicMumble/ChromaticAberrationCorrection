public class MyNum<T>{
    private T value;

    public MyNum(T value) {
        this.value = value;
    }

    public T get() {
        return value;
    }

    public void set(T anotherValue) {
        value = anotherValue;
    }

    @Override
    public String toString() {
        return value.toString();
    }

    @Override
    public boolean equals(Object obj) {
        return value.equals(obj);
    }
}
