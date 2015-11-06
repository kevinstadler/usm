(ns ^{:author "Kevin Stadler"}
  usm.core
  "Basic simulation functions for the utterance selection model of language change"
  (:refer-clojure :exclude [rand rand-int rand-nth])
  (:require
    [usm.measures :refer :all]
    [usm.random :refer :all]
    [clojure.java.io :refer [as-file delete-file writer]]
    [clojure.string :refer [join]]
    [incanter.stats :refer [mean]]))

(def pairsamplesperagent 100000)

; if every agent in the population comes this close to convergence on one of
; the variants we call the run off
(def epsilon 1.0E-6)

(defn limit [^double x]
  (min 1.0 (max 0.0 x)))

(defn ewma [^double alpha x y]
  (+ (* (- 1.0 alpha) x) (* alpha y)))

(defn apply-regularisation
  "Apply regularisation as in Blythe & Croft 2012 - parameter a should remain
  in [-2 1], positive a regularises, negative a 'de-regularises'"
  [^double a ^double x]
  (* x (inc (* a (- 1 x) (dec (* 2 x))))))

(defn apply-mrepl
  "Apply multiplicative replicator bias from Blythe & Croft 2012"
  [^double b ^double y]
  (case y
    (0.0 1.0) y
    (limit (* y (inc b)))))

(defn apply-repl [^double b ^double y]
  "Apply symmetric replicator bias (negative b is exact inverse)"
  (case y
    (0.0 1.0) y
    (limit (+ y b))))

(defn stochastic-prod
  "Return a random fraction B(n,p)/n drawn from a Binomial with n=T and p=x"
  [T ^double x]
  (/ (rand-binom T x) T))

(defn mutating-stochastic-prod
  "Return a random fraction B(n,p)/n drawn from a Binomial with n=T and p=x,
   with an additional possibility for tokens to mutate in production: m1 is the
   probability of producing an incoming variant when an old one was sampled, m2
   the probability of producing an old variant when an incoming one was sampled."
  [T ^double m1 ^double m2 ^double x]
  ; might also be possible to calculate modified p directly?
  (let [n (rand-binom T x)]
    (/ (+ n (rand-binom (- T n) m1) (- (rand-binom n m2))) T)))

(defn symmetric-mutating-stochastic-prod
  [T ^double m ^double x]
  (mutating-stochastic-prod T m m x))

(defn run-simulation
  [x0s pairsamplefun prodfun biasfun updatefun gammafun]
  (let [n (count x0s)
        population (make-array Double/TYPE 3 n) ; 3rd column for number of interactions that agent engaged in
        measures (make-array Double/TYPE (inc pairsamplesperagent) (* 3 (count (measurefn (repeat n 0)))))]
    ; init
    (dorun
      (map-indexed #(do (aset-double population 0 %1 %2)
                        (aset-double population 1 %1 %2))
                   x0s))
    ; momentum measures should all be zero, so store only initial x value measures
    (dorun (map-indexed #(aset-double measures 0 %1 %2) (measurefn (first population))))
    ; run
    (loop [i 1] ; count up to (at most) pairsamplesperagent
      (dotimes [j (/ n 2)] ; one interaction means two agent updates
        (let [[a1 a2] (pairsamplefun)
              x1 (aget population 0 a1)
              x2 (aget population 0 a2)
              g1 (aget population 1 a1)
              g2 (aget population 1 a2)
              y1 (prodfun x1)
              y2 (prodfun x2)
              fy1 (biasfun (- g2 x2) y1)
              fy2 (biasfun (- g1 x1) y2)]
          (aset-double population 0 a1 (updatefun x1 fy2))
          (aset-double population 0 a2 (updatefun x2 fy1))
          (aset-double population 1 a1 (gammafun g1 fy2))
          (aset-double population 1 a2 (gammafun g2 fy1))
          (aset-double population 2 a1 (inc (aget population 2 a1)))
          (aset-double population 2 a2 (inc (aget population 2 a2)))))
      (let [xs (first population)
            momentums (map - (second population) (first population))
            newmeasures (concat (measurefn xs) (measurefn momentums) (measurefn (nth population 2)))]
        (dorun (map-indexed #(aset-double measures i %1 %2) newmeasures))
        (if (or (== i pairsamplesperagent)
              (< (nth newmeasures 2) epsilon)
              (> (nth newmeasures 1) (- 1.0 epsilon)))
          (take (inc i) measures)
          (recur (inc i)))))))

(defn run-simulation-no-gamma
  [x0s pairsamplefun prodfun biasfun updatefun]
  (let [n (count x0s)
        population (make-array Double/TYPE 2 n) ; 2nd column for number of interactions that agent engaged in
        measures (make-array Double/TYPE (inc pairsamplesperagent) (* 2 (count (measurefn (repeat n 0)))))]
    (dorun (map-indexed #(aset-double population 0 %1 %2) x0s))
    (dorun (map-indexed #(aset-double measures 0 %1 %2) (measurefn (first population))))
    ; run
    (loop [i 1] ; count up to (at most) pairsamplesperagent
      (dotimes [j (/ n 2)] ; one interaction means two agent updates
        (let [[a1 a2] (pairsamplefun)
              x1 (aget population 0 a1)
              x2 (aget population 0 a2)
              y1 (prodfun x1)
              y2 (prodfun x2)
              fy1 (biasfun y1)
              fy2 (biasfun y2)]
          (aset-double population 0 a1 (updatefun x1 fy2))
          (aset-double population 0 a2 (updatefun x2 fy1))
          (aset-double population 1 a1 (inc (aget population 1 a1)))
          (aset-double population 1 a2 (inc (aget population 1 a2)))))
      (dorun (map-indexed #(aset-double measures i %1 %2) (concat (measurefn (first population)) (measurefn (second population)))))
      (if (or (== i pairsamplesperagent)
              (< (aget measures i 2) epsilon)
              (> (aget measures i 1) (- 1.0 epsilon)))
        (take (inc i) measures)
        (recur (inc i))))))

(defn write-conditions
  "Log all parameter combinations into a file readable by R"
  [prefix parameters]
  (.mkdir (as-file prefix))
  (let [f (as-file (str prefix "conditions.R"))]
    (when (not (.exists f))
      (with-open [w (writer f)]
        (doseq [[k v] parameters]
          (.write w (subs (str k " <- c(" (join ", " v) ")") 1))
          (.newLine w))))))

(defn write-data
  "Write simulation data to space-separated file (with column names header)"
  [filename data]
  (with-open [file (writer filename)]
    (.write file (join " " (mapcat #(map (fn [measure] (str measure %)) measurenames) observedvars)))
    (.newLine file)
    (doseq [line data] (.write file (join " " line)) (.newLine file))))
