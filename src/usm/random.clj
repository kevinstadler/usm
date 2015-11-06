(ns ^{:author "Kevin Stadler"}
  usm.random
  "Reproducible random sampling by thread-local random engine rebinding"
  (:refer-clojure :exclude [rand rand-int rand-nth])
  (:import
    (cern.jet.random.tdouble Binomial DoubleUniform)
    (cern.jet.random.tdouble.engine DoubleMersenneTwister)))

; the randomisation engine and distributions
(def ^:dynamic gen (DoubleMersenneTwister.))
(def ^:dynamic unif (DoubleUniform. gen))
(def ^:dynamic binom (Binomial. 1 0.5 gen))

(defmacro with-random-seed
  "Macro for temporarily resetting the current thread's Mersenne Twister
   randomisation engine with the given integer seed"
  [seed & body]
  ; bindings happen in parallel so do the engine first, distributions after
  `(binding [gen (DoubleMersenneTwister. ~seed)]
     (binding [unif (DoubleUniform. gen)
               binom (Binomial. 1 0.5 gen)]
     ~@body)))

(defn rand []
  (.nextDouble unif))

(defn rand-int
  "Uniformly sample an integer from [0...n-1]"
  [n]
  (.nextIntFromTo unif 0 (dec n)))

(defn rand-nth [coll]
  (nth coll (rand-int (count coll))))

(defn rand-int-outwith
  "Uniformly sample i integers from [0...n-1], excluding those given in the taboo
  list. This function will hang if i > n - |taboo|"
  ([n taboo i]
    (reduce (fn [rs _] (conj rs (rand-int-outwith n (concat taboo rs)))) [] (range i)))
  ([n taboo]
    (first (drop-while #(some (partial = %) taboo) (repeatedly #(.nextIntFromTo unif 0 (dec n)))))))
; sample 3 numbers from [0...4] without 1 and 3:
;(with-random-seed 1 (repeatedly 10 #(rand-int-outwith 5 [1 3] 3)))

(defn sample-pair
  "Uniformly sample two distinct integers from [0...n-1]"
  [n]
  (let [a1 (.nextIntFromTo unif 0 (dec n))
        a2 (mod (+ 1 a1 (.nextIntFromTo unif 0 (- n 2))) n)]
    [a1 a2]))

(defn permute [coll]
  (let [n (count coll)
        rs (conj (apply vector (map #(.nextIntFromTo unif 0 (dec %)) (range n 1 -1))) 0)
        is (reduce
             (fn [c idx] (let [[f r] (split-at idx c)] (concat f (map #(if (>= % (last f)) (inc %) %) r))))
             rs
             (range (dec n) 0 -1))]
    (map (partial nth coll) is)))
;(permute (range 5))
;(permute "abcde")
;(= (with-random-seed 1 (permute (range 5)))
;   (with-random-seed 1 (map #(- (int %) 97) (permute "abcde"))))

(defn rand-int-weighted
  "Sample a single integer from [0...(count weights)-1] using the given weights"
  [weights]
  (let [cumsum (reductions + weights)
        r (.nextDoubleFromTo unif 0.0 (last cumsum))]
    (first (keep-indexed (fn [i sm] (when (<= r sm) i)) cumsum))))

(defn rand-ints-weighted
  "Sample (up to) n distinct integers from [0...(count weights)-1] using the given weights"
  [weights n]
  (if (<= (count (keep pos? weights)) n)
    (permute (keep-indexed #(when (pos? %2) %1) weights))
    (let [cumsum (reductions + weights)]
      (doall ; force immediate evaluation to stop lazy evaluation bypassing random engine rebinding
        (take n
          (distinct
            (repeatedly
              #(let [r (.nextDoubleFromTo unif 0.0 (last cumsum))]
                (first (keep-indexed (fn [i sm] (when (<= r sm) i)) cumsum))))))))))
;(= [3 4] (sort (rand-ints-weighted [0 0 0 1 1] 2)))
;(= [0 4] (sort (rand-ints-weighted [1 0 0 0 1] 2)))

(defn rand-nth-weighted [weights coll]
  (nth coll (rand-int-weighted weights)))

(defn rand-binom
  "Sample a Binomial from Parallelcolt, handling special cases of n and p."
  [n ^double p]
  (if (zero? n)
    0
    (case p
      0.0 0
      1.0 n
      (.nextInt binom n p))))
