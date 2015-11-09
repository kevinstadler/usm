(ns ^{:author "Kevin Stadler"}
  usm.experiments.momentum
  "Run batches of the momentum-based model of language change"
  (:gen-class)
  (:refer-clojure :exclude [rand rand-int rand-nth])
  (:require
    [usm.core :refer :all]
    [usm.momentum :refer :all]
    [usm.random :refer [sample-pair]]
    [clojure.java.io :refer [as-file delete-file writer]]))

(def conditions {
  :a [0.01]
  :t [2 3 4 5]
  :x0 [0.01 0.02 0.03]
  :i [0 1]
  :b [0 0.5 1 1.5 2 2.5]
  :m [1.5 2 2.5 3 3.5 4]
  :N [2 5 10 20 30 50]; 100 200 300])
})

(defn run-conditions [prefix]
  (doseq [alpha (:a conditions)
          T (:t conditions)
          x0 (:x0 conditions)
          mult (:m conditions)
          inverted (:i conditions)
          b (:b conditions)
          N (:N conditions)]
    (let [filename (str prefix "N" N "b" b "T" T "x0" x0 "i" inverted "a" alpha "m" mult)]
      (run-continuous-momentum (partial sample-pair N) b T x0 inverted N alpha mult filename))))

; to cleanly abort the simulations after the next batch of runs has finished,
; simply create a file named 'stop' in the directory where simulations are run
(def stop-file (as-file "stop"))
(def parallelruns 8)

(defn -main [& args]
  (let [timestamp (if (> (count args) 2) (Long/parseLong (nth args 2)) (System/currentTimeMillis))
        prefix (str "results-clj-" timestamp "/")
        start (Integer/parseInt (first args))
        end (Integer/parseInt (or (second args) (str (+ start parallelruns))))]
    (println "Using" (if (> (count args) 2) "provided" "current") "timestamp" timestamp "as seed for random number generation")
    (loop [i start]
      (when (and (<= i end) (not (.exists stop-file)))
        (println "Running" (apply * (map count (vals conditions))) "conditions" parallelruns "times each")
        (write-conditions prefix conditions)
        (let [run-prefixes (map #(str prefix % "/") (range i (min (inc end) (+ i parallelruns))))
              status-file (str "running-" i "-" (dec (+ i parallelruns)))]
          ; create target dirs
          (doseq [pref run-prefixes]
            (when (not (.mkdir (as-file pref)))
              (println "Output directory" pref "already exists, aborting!")
              (System/exit 1)))
          (spit status-file "")
          (time (dorun (pmap run-conditions run-prefixes)))
          (delete-file status-file))
        (recur (+ i parallelruns))))))

; lein uberjar
; nohup time java -jar usm-1.0.0-standalone.jar 1 8 &
