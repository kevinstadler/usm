(ns ^{:author "Kevin Stadler"}
  usm.momentum
  "Basic simulation functions for the momentum-based model of language change"
  (:refer-clojure :exclude [rand rand-int rand-nth])
  (:require
    [usm.core :refer :all]
    [usm.measures :refer :all]
    [usm.random :refer [sample-pair with-random-seed]]
    [clojure.java.io :refer [as-file delete-file writer]]))

(defn apply-momentum
  "Apply symmetric version of Gureckis & Goldstone momentum bias"
  [^double b m ^double y]
  (case y
    (0.0 1.0) y
    (limit (+ y (* b m)))))

(defn apply-updown-momentum
  "Apply simple binary up/down momentum bias"
  [^double b m ^double y]
  (case y
    (0.0 1.0) y
    (limit (+ y (* b (Math/signum m))))))

(defn stepstomaxdifference [alpha gamma]
  (/ (Math/log (/ alpha gamma)) (- alpha gamma)))

(defn expdecay [lambda t]
  (Math/exp (- (* lambda t))))

(defn maxdecaydifference [alpha gamma]
  (let [t (stepstomaxdifference alpha gamma)]
    (- (expdecay alpha t) (expdecay gamma t))))

(defn normaliseb [alpha gamma b]
  (/ b (maxdecaydifference alpha gamma)))
;(normaliseb 0.01 0.015 0.1)

(defn run-momentum [agent-samplefun biasfun T x0 inverted N alpha mult filename]
  (with-random-seed (.hashCode filename)
    (let [sample-data (partial stochastic-prod T)
          updatefun (partial ewma alpha)
          gammafun (partial ewma (* alpha mult))
          x0fun (constantly (if (zero? inverted) x0 (- 1 x0)))]
      (write-data filename (run-simulation (repeatedly N x0fun) agent-samplefun sample-data biasfun updatefun gammafun)))))

(defn run-continuous-momentum [agent-samplefun b T x0 inverted N alpha mult filename]
  (let [normalisedb (normaliseb alpha (* alpha mult) b)]
    (run-momentum agent-samplefun (partial apply-momentum normalisedb) T x0 inverted N alpha mult filename)))

; binary up/down momentum: bias strength directly comparable to the additive replicator selection one
(defn run-simple-momentum [agent-samplefun b T x0 inverted N alpha mult filename]
  (with-discrete-momentum-measures
    (run-momentum agent-samplefun (partial apply-updown-momentum b) T x0 inverted N alpha mult filename)))

;(run-simple-momentum #(sample-pair 20) 0.005 2 0.01 0 20 0.01 2.5 "foo2.txt")
