from tensorflow.keras import models, layers, utils, backend as K
import NeuralNetworks.VisualizationNN as vn
import matplotlib.pyplot as plt
import numpy as np


def TrainANN(Architecture, InputsSamples, OutputSamples, NormalizationOption, ActivationFunc, DropoutOption, batchSize = 320, numEpochs = 3000):

    # normalization 
    normalizer = layers.Normalization(axis=-1)
    normalizer.adapt(np.array(InputsSamples))
    print(normalizer.mean.numpy())
    myLayers = [normalizer]

    # Hidden layers
    numHidden = Architecture.size
    for i in range(numHidden):
        hiddenLayerName = "Hidden_{index}".format(index = i)
        hiddenLayer = layers.Dense(name = hiddenLayerName, units = Architecture[i], activation = ActivationFunc)
        myLayers.append(hiddenLayer)

        # Dropout layer to avoid overfitting
        if DropoutOption == 1:
            dropOutLayerName = "Dropout_{index}".format(index = i)
            dropOutLayer = layers.Dropout(name = dropOutLayerName, rate = 0.05)
            myLayers.append(dropOutLayer)

    # Output layer
    shapeOutSam = OutputSamples.shape
    if len(shapeOutSam) == 1:
        numOutputs = 1
    else:
        numOutputs = OutputSamples.shape[1]
    outputLayer = layers.Dense(name = "output", units = numOutputs, activation='linear')
    myLayers.append(outputLayer)

    # n_features =30
    # hidden_1 = layers.Dense(name="h1", units=int(round((n_features+1)/2)), activation='relu')
    # drop_1 = layers.Dropout(name="drop1", rate=0.2)
    # hidden_2 = layers.Dense(name="h2", units=int(round((n_features+1)/2)), activation='relu')
    # drop_2 = layers.Dropout(name="drop2", rate=0.2)
    # out = layers.Dense(name="output", units=1, activation='sigmoid')

    # myLayers2 = [hidden_1, drop_1, hidden_2, drop_2, out]




    # Create ANN
    myANN = models.Sequential(name = "DeepNN", layers = myLayers)
    myANN.summary()
    # vn.visualize_nn(myANN, description=True, figsize=(10,8))


    # compile the neural network
    myANN.compile(optimizer = "Adam", loss = "mse", metrics = ["mae"])
    # myANN.compile(optimizer = "Adadelta", loss = "mse", metrics = ["mse"])

    # Normalize the output samples as well
    OutMax = np.amax(OutputSamples, axis=0)
    OutMin = np.amin(OutputSamples, axis=0)
    OutDelta = OutMax - OutMin

    if NormalizationOption == 0:
        OutputSamples_Train = OutputSamples
    elif NormalizationOption == 1:
        OutputSamples_Norm = (OutputSamples - OutMin) / OutDelta  # Normalize to (0, 1)
        OutputSamples_Train = OutputSamples_Norm
    elif NormalizationOption == 2:
        OutputSamples_Norm = (OutputSamples - OutMin) / OutDelta * 2 - 1 # Normalize to (-1, 1)
        OutputSamples_Train = OutputSamples_Norm
    # train/validation
    training = myANN.fit(x = InputsSamples, y = OutputSamples_Train, batch_size=batchSize, epochs=numEpochs, shuffle=True, verbose=0, validation_split=0.25)


    # plot Training
    metrics = [k for k in training.history.keys() if ("loss" not in k) and ("val" not in k)]    
    fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(15,3))
            
    ax[0].set(title="Training")    
    ax11 = ax[0].twinx()    
    ax[0].plot(training.history['loss'], color='black')    
    ax[0].set_xlabel('Epochs')    
    ax[0].set_ylabel('Loss', color='black')    
    for metric in metrics:        
        ax11.plot(training.history[metric], label=metric)    
    ax11.set_ylabel("Score", color='steelblue')    
    ax11.legend()
            
    ## validation    
    ax[1].set(title="Validation")    
    ax22 = ax[1].twinx()    
    ax[1].plot(training.history['val_loss'], color='black')    
    ax[1].set_xlabel('Epochs')    
    ax[1].set_ylabel('Loss', color='black')    
    for metric in metrics:          
        ax22.plot(training.history['val_'+metric], label=metric)    
    ax22.set_ylabel("Score", color="steelblue")    
    plt.show(block=False)

    return (myANN, OutMax, OutMin, OutDelta)

