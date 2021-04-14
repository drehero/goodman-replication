# Konzeption Paper

- Zieljournal: The American Statistician
- andere Möglichkeiten: Rescience C, ganz evtl. Plos One


## Introduction
- Blabla p-values und statistische Inferenz werden kritisiert, American Statistician mit special issue zu der Thematik
- In einem Paper vergleicht Goodman verschiedene (frequentistische Testverfahren), er findet, dass…
- wir replizieren die Ergebnisse von Goodman in einem engen Sinne, stellen Software Code zur Verfügung (github) und ergänzen die Ergebnisse um weitere statistische Testverfahren
- "We believe that an openly available code repository replicating the results of Goodman’s paper can be helpful to the scientifi
community. We therefore implemented the model and analysis scripts using R..."

## Methods
- was macht Goodman? Statistische Tests von ihm und zusätzliche von uns sowie Simulationssetting beschreiben (auch mögliche Äquivalenz beschreiben, wie z.B. dass die Methode von Betensky gleich ist zur Interval-based method); Erweiterungen des Simulationssetting wie in Goodman’s Appendix beschrieben bzw. vorgeschlagen könnten hier oder in einem Appendix eingefügt werden)

## Welche Methoden fügen wir hinzu?
- „flat prior“ um 5%-Fehlerrate unter H0 einzuhalten wie von Arne vorgeschlagen
- ansonsten aus dem special issue des American Statistician:
-- Equivalence Testing
-- Blume’s method 
-- Evtl. Für p-werte und equivalence testing Werte zwischen 0.5 und 5% als suggestive bezeichnen (also quasi rausschmeißen)
-- Methoden von Colquhoun (False positive risk) und Matthews (Intrinsic Analysis of Credibility) hinzufügen?
-- ...weitere?

## Results
- Replikation erfolgreich, Ergebnisse beschreiben
- Ergebnisse zusätzlicher Methoden beschreiben

## Conclusions

